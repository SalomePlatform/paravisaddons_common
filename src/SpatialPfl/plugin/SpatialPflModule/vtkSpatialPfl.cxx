// Copyright (C) 2021  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
// Author : Anthony Geay (EDF R&D)

#include "vtkSpatialPfl.h"

#include <vtkAdjacentVertexIterator.h>
#include <vtkAlgorithmOutput.h>
#include <vtkCell.h>
#include <vtkCellType.h>
#include <vtkCharArray.h>
#include <vtkDataArraySelection.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationDataObjectKey.h>
#include <vtkInformationStringKey.h>
#include <vtkInformationVector.h>
#include <vtkLineSource.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkResampleWithDataSet.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTable.h>
#include <vtkTimeStamp.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVariantArray.h>
#include "vtkPolyData.h"

#include <deque>
#include <map>
#include <sstream>

vtkStandardNewMacro(vtkSpatialPfl);

///////////////////

template <class T>
class AutoPtr
{
public:
  AutoPtr(T* ptr = nullptr)
    : _ptr(ptr)
  {
  }
  ~AutoPtr() { destroyPtr(); }
  AutoPtr& operator=(T* ptr)
  {
    if (_ptr != ptr)
    {
      destroyPtr();
      _ptr = ptr;
    }
    return *this;
  }
  T* operator->() { return _ptr; }
  const T* operator->() const { return _ptr; }
  T& operator*() { return *_ptr; }
  const T& operator*() const { return *_ptr; }
  operator T*() { return _ptr; }
  operator const T*() const { return _ptr; }

private:
  void destroyPtr() { delete[] _ptr; }

private:
  T* _ptr;
};

class MZCException : public std::exception
{
public:
  MZCException(const std::string& s)
    : _reason(s)
  {
  }
  virtual const char* what() const throw() { return _reason.c_str(); }
  virtual ~MZCException() throw() {}

private:
  std::string _reason;
};

void ExtractInfo(vtkInformationVector* inputVector, vtkUnstructuredGrid*& usgIn)
{
  vtkInformation* inputInfo(inputVector->GetInformationObject(0));
  vtkDataSet* input = nullptr;
  vtkDataSet* input0 = vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkMultiBlockDataSet* input1 = vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT()));
  if (input0)
  {
    input = input0;
  }
  else
  {
    if (!input1)
    {
      throw MZCException(
        "Input dataSet must be a DataSet or single elt multi block dataset expected !");
    }
    if (input1->GetNumberOfBlocks() != 1)
    {
      throw MZCException("Input dataSet is a multiblock dataset with not exactly one block ! Use "
                         "MergeBlocks or ExtractBlocks filter before calling this filter !");
    }
    vtkDataObject* input2 = input1->GetBlock(0);
    if (!input2)
    {
      throw MZCException("Input dataSet is a multiblock dataset with exactly one block but this "
                         "single element is NULL !");
    }
    vtkDataSet* input2c = vtkDataSet::SafeDownCast(input2);
    if (!input2c)
    {
      throw MZCException(
        "Input dataSet is a multiblock dataset with exactly one block but this single element is "
        "not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    }
    input = input2c;
  }
  if (!input)
  {
    throw MZCException("Input data set is NULL !");
  }
  usgIn = vtkUnstructuredGrid::SafeDownCast(input);
  if (!usgIn)
  {
    throw MZCException("Input data set is not an unstructured mesh ! This filter works only on "
                       "unstructured meshes !");
  }
}

////////////////////

vtkSpatialPfl::vtkSpatialPfl()
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int vtkSpatialPfl::RequestInformation(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // std::cerr << "########################################## vtkSpatialPfl::RequestInformation
  // ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid* usgIn(0);
    ExtractInfo(inputVector[0], usgIn);
    vtkInformation* info(outputVector->GetInformationObject(0));
  }
  catch (MZCException& e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkSpatialPfl::RequestInformation : " << e.what()
        << std::endl;
    if (this->HasObserver("ErrorEvent"))
    {
      this->InvokeEvent("ErrorEvent", const_cast<char*>(oss.str().c_str()));
    }
    else
    {
      vtkOutputWindowDisplayErrorText(const_cast<char*>(oss.str().c_str()));
    }
    vtkObject::BreakOnError();
    return 0;
  }
  return 1;
}

void buildTableFrom(vtkPolyData* table, const std::vector<std::vector<double> >& valuesByColumn, const std::vector<std::string>& columnNames)
{
  vtkPointData *pd(table->GetPointData());
  std::size_t sz(valuesByColumn.size());
  if (sz != columnNames.size())
  {
    throw MZCException("Sizes of vectors mismatches !");
  }
  if (sz == 0)
  {
    return;
  }
  std::size_t nbSamples(valuesByColumn[0].size());
  for (int i = 0; i < sz; i++)
  {
    vtkSmartPointer<vtkDoubleArray> arr(vtkSmartPointer<vtkDoubleArray>::New());
    arr->SetName(columnNames[i].c_str());
    if (nbSamples != valuesByColumn[i].size())
    {
      std::ostringstream oss;
      oss << "Sizes of vectors " << i << " mismatches with size of others " << nbSamples << " !";
      throw MZCException(oss.str());
    }
    arr->SetNumberOfTuples(nbSamples);
    arr->SetNumberOfComponents(1);
    double* pt(arr->GetPointer(0));
    std::copy(valuesByColumn[i].begin(), valuesByColumn[i].end(), pt);
    pd->AddArray(arr);
  }
}
template<class T>
struct VTKArrayTraits
{
};

template<>
struct VTKArrayTraits<double>
{
  using Type = vtkDoubleArray;
};

template<>
struct VTKArrayTraits<float>
{
  using Type = vtkFloatArray;
};

template<class T>
void FromVTKArrayComputeCurvAbsInternal(typename VTKArrayTraits<T>::Type *data, std::vector<double>& ret, std::vector<double>& Xcoords, std::vector<double>& Ycoords)
{
  vtkIdType nbTuples(data->GetNumberOfTuples()), nbComp(data->GetNumberOfComponents());
  ret.resize(nbTuples); Xcoords.resize(nbTuples); Ycoords.resize(nbTuples);
  const T* pt(data->GetPointer(0));
  ret[0] = 0.; Xcoords[0] = pt[0]; Ycoords[0] = pt[1];
  for (vtkIdType i = 1; i < nbTuples; i++)
  {
    double val(0.);
    for (vtkIdType j = 0; j < nbComp; j++)
    {
      double delta(pt[nbComp * (i - 1) + j] - pt[nbComp * i + j]);
      val += delta * delta;
    }
    Xcoords[i] = pt[nbComp * i + 0]; Ycoords[i] = pt[nbComp * i + 1];
    ret[i] = ret[i - 1] + sqrt(val);
  }
}

std::vector<double> FromVTKArrayComputeCurvAbs(vtkDataArray* data, std::vector<double>& Xcoords, std::vector<double>& Ycoords)
{
  vtkIdType nbTuples(data->GetNumberOfTuples()), nbComp(data->GetNumberOfComponents());
  if (nbTuples < 1)
  {
    throw MZCException("FromVTKArrayComputeCurvAbs : internal error 1 !");
  }
  std::vector<double> ret;
  vtkDoubleArray* d1(vtkDoubleArray::SafeDownCast(data));
  if (d1)
  {
    FromVTKArrayComputeCurvAbsInternal<double>(d1,ret,Xcoords,Ycoords);
    return ret;
  }
  vtkFloatArray* d2(vtkFloatArray::SafeDownCast(data));
  if (d2)
  {
    FromVTKArrayComputeCurvAbsInternal<float>(d2,ret,Xcoords,Ycoords);
    return ret;
  }
  throw MZCException("FromVTKArrayComputeCurvAbs : internal error 2 !");
}

static std::string GetReprDependingPos(const std::string& origName, int blockId, int nbBlocks)
{
  if( nbBlocks == 1 )
    return origName;
  std::ostringstream oss;
  oss << origName << "_" << blockId;
  return oss.str();
}

void FillPolyDataInstance(vtkPolyData *ds, const std::vector<double>& xs, const std::vector<double>& ys)
{
  vtkNew<vtkDoubleArray> coords;
  std::size_t nbPts(xs.size());
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(nbPts);
  double *coordsPt(coords->GetPointer(0));
  vtkNew<vtkPoints> pts;
  pts->SetData(coords);
  ds->SetPoints(pts);
  //
  vtkNew<vtkCellArray> cc;
  cc->AllocateExact(1,nbPts);
  std::unique_ptr<vtkIdType[]> conn(new vtkIdType[nbPts]);
  //
  for(std::size_t iPt = 0 ; iPt < nbPts ; ++iPt)
    {
      coordsPt[3*iPt] = xs[iPt] ; coordsPt[3*iPt+1] = ys[iPt] ; coordsPt[3*iPt+2] = 0.;
      conn[iPt] = iPt;
    }
  //
  cc->InsertNextCell(nbPts,conn.get());
  ds->SetLines(cc);
}

void buildTableFromPolyData(vtkMultiBlockDataSet* table, vtkPolyData* ds, int blockId, int nbBlocks)
{
  vtkNew<vtkPolyData> eltInTable;
  vtkPoints* pts(ds->GetPoints());
  if (!pts)
  {
    throw MZCException("buildTableFromPolyData : internal error 2 !");
  }
  vtkDataArray* data(pts->GetData());
  if (!data)
  {
    throw MZCException("buildTableFromPolyData : internal error 3 !");
  }
  //
  std::vector<double> Xcoords, Ycoords;
  std::vector<double> xs(FromVTKArrayComputeCurvAbs(data, Xcoords, Ycoords));
  //
  FillPolyDataInstance(eltInTable,Xcoords,Ycoords);
  //
  vtkCellArray *cd(ds->GetVerts()), *cc(ds->GetLines());
  //
  vtkDataSetAttributes* dsa(ds->GetPointData());
  if (!dsa)
  {
    throw MZCException("buildTableFromPolyData : no point data !");
  }
  int nba = dsa->GetNumberOfArrays();
  //
  std::vector<std::vector<double> > valuesByColumn(1);
  std::vector<std::string> columnNames(1);
  valuesByColumn[0] = xs;
  columnNames[0] = "Curv Abscissa";
  //
  for (int i = 0; i < nba; i++)
  {
    vtkDataArray* arr(dsa->GetArray(i));
    std::vector<double> tmp(arr->GetNumberOfTuples());
    if (arr->GetNumberOfComponents() != 1)
    {
      continue;
    }
    std::string name(GetReprDependingPos(arr->GetName(),blockId,nbBlocks));
    vtkDoubleArray* arr1(vtkDoubleArray::SafeDownCast(arr));
    if (!arr1)
    {
      vtkFloatArray* arr2(vtkFloatArray::SafeDownCast(arr));
      if (!arr2)
      {
        continue;
      }
      const float* pt(arr2->GetPointer(0));
      std::copy(pt, pt + arr->GetNumberOfTuples(), tmp.begin());
    }
    else
    {
      const double* pt(arr1->GetPointer(0));
      std::copy(pt, pt + arr1->GetNumberOfTuples(), tmp.begin());
    }
    valuesByColumn.push_back(tmp);
    columnNames.push_back(name);
  }
  // EDF21757 - Ajout de X et Y
  valuesByColumn.push_back(Xcoords); columnNames.push_back("X");
  valuesByColumn.push_back(Ycoords); columnNames.push_back("Y");
  //
  buildTableFrom(eltInTable, valuesByColumn, columnNames);
  table->SetBlock(blockId,eltInTable);
}

vtkPolyData *ExtractTo(vtkDataObject *source)
{
  vtkMultiBlockDataSet *sourceMB(vtkMultiBlockDataSet::SafeDownCast(source));
  if(sourceMB)
  {
    if(sourceMB->GetNumberOfBlocks() != 1)
    {
      std::ostringstream oss; oss << "Internal error ! Number of blocks of MultiBlockDataSet source must be equal to 1 ! Here : " << sourceMB->GetNumberOfBlocks();
      throw MZCException(oss.str());
    }
    vtkDataObject *source20(sourceMB->GetBlock(0));
    vtkPolyData *source20c(vtkPolyData::SafeDownCast(source20));
    if(!source20c)
    {
      throw MZCException("Internal error ! source is a mono block MultiBlockDataSet but this block is not a vtkDataSet !");
    }
    return source20c;
  }
  vtkPolyData* source2(vtkPolyData::SafeDownCast(source));
  return source2;
}

static int GetNumberOfBlocs(vtkDataObject *ds)
{
  if(!ds)
    throw MZCException("vtkSedimentDeposit  SplitSingleMultiBloc : nullptr !");
  vtkMultiBlockDataSet *ds0(vtkMultiBlockDataSet::SafeDownCast(ds));
  if(!ds0)
  {
    vtkPolyData *ds00(vtkPolyData::SafeDownCast(ds));
    if(!ds00)
      throw MZCException("vtkSedimentDeposit  SplitSingleMultiBloc : neither a vtkMultiBlockDataSet nor a vtkPolyData !");
    return 1;
  }
  return ds0->GetNumberOfBlocks();
}

static vtkPolyData *SplitSingleMultiBloc(vtkDataObject *res, int blockId)
{
  vtkPolyData* res2(vtkPolyData::SafeDownCast(res));
  if (!res2)
  {
    vtkMultiBlockDataSet *res2c(vtkMultiBlockDataSet::SafeDownCast(res));
    if(!res2c)
    {
      throw MZCException("Internal error ! unexpected returned of resample filter !");
    }
    if(blockId >= res2c->GetNumberOfBlocks())
    {
      std::ostringstream oss; oss << "Internal error ! Number of blocks of MultiBlockDataSet must be equal < " << blockId << " ! Here : " << res2c->GetNumberOfBlocks();
      throw MZCException(oss.str());
    }
    vtkDataObject *res20(res2c->GetBlock(blockId));
    vtkPolyData *res20c(vtkPolyData::SafeDownCast(res20));
    if(!res20c)
    {
      throw MZCException("Internal error ! resample filter returned a mono block MultiBlockDataSet but this block is not a vtkPolyData !");
    }
    return res20c;
  }
  if(blockId!=0)
    throw MZCException("Internal error ! 0 expected !");
  return res2;
}

int vtkSpatialPfl::RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // std::cerr << "########################################## vtkSpatialPfl::RequestData
  // ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid* usgIn(0);
    ExtractInfo(inputVector[0], usgIn);
    //
    vtkInformation* sourceInfo(inputVector[1]->GetInformationObject(0));
    vtkDataObject* source(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkDataSet* source2(vtkDataSet::SafeDownCast(source));

    int nbBlocks2(GetNumberOfBlocs(source));
    vtkNew<vtkMultiBlockDataSet> sampleLineMB;
    sampleLineMB->SetNumberOfBlocks(nbBlocks2);

    //vtkSmartPointer<vtkDataSet> sampleLine;
    if (this->ResampleInput)
    {
      for(int iBlock = 0 ; iBlock < nbBlocks2 ; ++iBlock)
      {
        vtkPolyData* pl = SplitSingleMultiBloc(source,iBlock);
        if (!pl)
        {
          vtkErrorMacro("The second input of this filter must be of type vtkPolyData.");
          return 1;
        }
        vtkSmartPointer<vtkDataSet> sampleLine;
        sampleLine.TakeReference(this->ResamplePolyLine(pl));
        sampleLineMB->SetBlock(iBlock,sampleLine);
      }
    }
    else
    {
      for(int iBlock = 0 ; iBlock < nbBlocks2 ; ++iBlock)
      {
        //sampleLine = ExtractTo(source);
        sampleLineMB->SetBlock(iBlock,SplitSingleMultiBloc(source,iBlock));
      }
    }
    //
    vtkNew<vtkResampleWithDataSet> probeFilter;
    probeFilter->SetInputData(sampleLineMB);
    probeFilter->SetSourceData(usgIn);
    probeFilter->Update();
    vtkDataObject* res(probeFilter->GetOutput());
    //
    int nbBlocks(GetNumberOfBlocs(res));
    vtkNew<vtkTable> table;
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkMultiBlockDataSet *output ( vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())) );
    output->SetNumberOfBlocks(nbBlocks);
    for(int blockId = 0 ; blockId < nbBlocks ; ++blockId)
    {
      vtkPolyData *res2(SplitSingleMultiBloc(res,blockId));
      buildTableFromPolyData(output, res2, blockId, nbBlocks);
    }
  }
  catch (MZCException& e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkSpatialPfl::RequestData : " << e.what() << std::endl;
    if (this->HasObserver("ErrorEvent"))
    {
      this->InvokeEvent("ErrorEvent", const_cast<char*>(oss.str().c_str()));
    }
    else
    {
      vtkOutputWindowDisplayErrorText(const_cast<char*>(oss.str().c_str()));
    }
    vtkObject::BreakOnError();
    return 0;
  }
  return 1;
}

void vtkSpatialPfl::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkSpatialPfl::SetSourceData(vtkDataObject* input)
{
  this->SetInputData(1, input);
}

void vtkSpatialPfl::SetSourceConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}

int vtkSpatialPfl::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
}

vtkPolyData* vtkSpatialPfl::ResamplePolyLine(vtkPolyData* pl)
{
  vtkPolyData* res = vtkPolyData::New();

  vtkNew<vtkLineSource> subdivision;
  subdivision->SetPoints(pl->GetPoints());
  subdivision->SetResolution(this->NumberOfSamples);
  subdivision->Update();
  res->ShallowCopy(subdivision->GetOutputDataObject(0));

  return res;
}

/*<Hints>
        <!-- View can be used to specify the preferred view for the proxy -->
        <View type="QuartileChartView" />
        </Hints>*/

//   /opt/cmake/3.6.2/bin/cmake -DCMAKE_INSTALL_PREFIX=/home/H87074/TMP117_HYDRAU/SpatialPfl_install
//   -DCONFIGURATION_ROOT_DIR=/opt/salome-conf/8.3.0 -DCMAKE_BUILD_TYPE=Debug
//   -DPYTHON_INCLUDE_DIR=/usr/include/python2.7
//   -DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython2.7.so ../SpatialPfl
/*
export LD_LIBRARY_PATH=/home/H87074/TMP117_HYDRAU/SpatialPfl_install/lib:${LD_LIBRARY_PATH}
export PV_PLUGIN_PATH=/home/H87074/TMP117_HYDRAU/SpatialPfl_install/lib:${PV_PLUGIN_PATH}
*/
