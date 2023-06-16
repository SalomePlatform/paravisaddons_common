// Copyright (C) 2021-2023  CEA/DEN, EDF R&D
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

#include "vtkDepthVsTime.h"

#include <vtkAdjacentVertexIterator.h>
#include <vtkAlgorithmOutput.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCharArray.h>
#include <vtkCompositeDataProbeFilter.h>
#include <vtkDataArraySelection.h>
#include <vtkDataObjectTreeIterator.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkFloatArray.h>
#include <vtkInEdgeIterator.h>
#include <vtkInformation.h>
#include <vtkInformationDataObjectKey.h>
#include <vtkInformationStringKey.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>
#include <vtkResampleWithDataSet.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkTable.h>
#include <vtkTimeStamp.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVariantArray.h>
#include <vtkWarpScalar.h>

#include <deque>
#include <map>
#include <set>
#include <sstream>

vtkStandardNewMacro(vtkDepthVsTime);

constexpr int NB_OF_DISCR_TO_DEDUCE_START_STOP = 1001;

///////////////////

template <class T>
class AutoPtr
{
public:
  AutoPtr(T *ptr = 0) : _ptr(ptr) {}
  ~AutoPtr() { destroyPtr(); }
  AutoPtr &operator=(T *ptr)
  {
    if (_ptr != ptr)
    {
      destroyPtr();
      _ptr = ptr;
    }
    return *this;
  }
  T *operator->() { return _ptr; }
  const T *operator->() const { return _ptr; }
  T &operator*() { return *_ptr; }
  const T &operator*() const { return *_ptr; }
  operator T *() { return _ptr; }
  operator const T *() const { return _ptr; }

private:
  void destroyPtr() { delete[] _ptr; }

private:
  T *_ptr;
};

class MZCException : public std::exception
{
public:
  MZCException(const std::string &s) : _reason(s) {}
  virtual const char *what() const throw() { return _reason.c_str(); }
  virtual ~MZCException() throw() {}

private:
  std::string _reason;
};

class vtkDepthVsTime::vtkInternal
{
public:
  vtkInternal() : _isInit(true) { _nbItems = std::numeric_limits<std::size_t>::max(); }
  void operate(double timeStep, vtkUnstructuredGrid *usgIn);
  void pushData(double timeStep, vtkPolyData *data);
  void fillGrid(vtkRectilinearGrid *grid) const;
  void scanCoordsOfDS(vtkUnstructuredGrid *usg, vtkPolyData *zePoint, int nbOfExpectedDiscr);
  static std::size_t CheckPts(vtkPointSet *usg, vtkDataArray *&arr);

private:
  void pushDataInit(double timeStep, vtkDataSetAttributes *dsa);
  void pushDataStd(double timeStep, vtkDataSetAttributes *dsa);

private:
  bool _isInit;
  std::size_t _nbItems;
  std::vector<std::string> _columnNames;
  std::vector<double> _time;
  // sizeof(_data)==sizeof(_columnNames)
  std::vector<std::vector<double>> _data;
  vtkSmartPointer<vtkPolyData> _ladder;
};

void ExtractInfo(vtkInformationVector *inputVector, vtkUnstructuredGrid *&usgIn)
{
  vtkInformation *inputInfo(inputVector->GetInformationObject(0));
  vtkDataSet *input(0);
  vtkDataSet *input0(vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  vtkMultiBlockDataSet *input1(vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  if (input0)
    input = input0;
  else
  {
    if (!input1)
      throw MZCException("Input dataSet must be a DataSet or single elt multi block dataset expected !");
    if (input1->GetNumberOfBlocks() != 1)
      throw MZCException("Input dataSet is a multiblock dataset with not exactly one block ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    vtkDataObject *input2(input1->GetBlock(0));
    if (!input2)
      throw MZCException("Input dataSet is a multiblock dataset with exactly one block but this single element is NULL !");
    vtkDataSet *input2c(vtkDataSet::SafeDownCast(input2));
    if (!input2c)
      throw MZCException("Input dataSet is a multiblock dataset with exactly one block but this single element is not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    input = input2c;
  }
  if (!input)
    throw MZCException("Input data set is NULL !");
  usgIn = vtkUnstructuredGrid::SafeDownCast(input);
  if (!usgIn)
    throw MZCException("Input data set is not an unstructured mesh ! This filter works only on unstructured meshes !");
}

////////////////////

vtkDepthVsTime::vtkDepthVsTime() : NbDiscrPtsAlongZ(10), NumberOfTimeSteps(0), IsExecuting(false), CurrentTimeIndex(0), Internal(NULL)
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

vtkDepthVsTime::~vtkDepthVsTime()
{
  delete this->Internal;
  this->Internal = NULL;
}

int vtkDepthVsTime::RequestInformation(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkDepthVsTime::RequestInformation ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid *usgIn(0);
    ExtractInfo(inputVector[0], usgIn);
    vtkInformation *inInfo(inputVector[0]->GetInformationObject(0));
    if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
    {
      this->NumberOfTimeSteps = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    }
    else
    {
      this->NumberOfTimeSteps = 0;
    }
    // The output of this filter does not contain a specific time, rather
    // it contains a collection of time steps. Also, this filter does not
    // respond to time requests. Therefore, we remove all time information
    // from the output.
    vtkInformation *outInfo(outputVector->GetInformationObject(0));
    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
    {
      outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    }
    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_RANGE()))
    {
      outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());
    }
    return 1;
  }
  catch (MZCException &e)
  {
    vtkErrorMacro(<< "Exception has been thrown in vtkDepthVsTime::RequestInformation : " << e.what());
    return 0;
  }
  return 1;
}

int vtkDepthVsTime::RequestUpdateExtent(vtkInformation *, vtkInformationVector **inputVector, vtkInformationVector *vtkNotUsed(outputVector))
{
  // vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkInformation *inInfo1 = inputVector[0]->GetInformationObject(0);

  // get the requested update extent
  double *inTimes = inInfo1->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  if (inTimes)
  {
    double timeReq = inTimes[this->CurrentTimeIndex];
    inInfo1->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), timeReq);
  }

  return 1;
}

std::string buildNameOfEntryFrom(const std::string &name, std::size_t id, std::size_t nbOfElt)
{
  if (nbOfElt == 0)
    throw MZCException("buildNameOfEntryFrom : nbElt == 0 !");
  if (nbOfElt == 1)
    return name;
  std::ostringstream oss;
  oss << name << "_" << id;
  return oss.str();
}

void fillFieldInGrid(vtkDataSet *ds, const std::vector<std::vector<double>> &valuesByColumn, std::size_t szX, std::size_t szY, const std::vector<std::string> &columnNames)
{
  vtkPointData *pd(ds->GetPointData());
  if (!pd)
    throw MZCException("fillFieldInGrid : internal error 1 !");
  std::size_t nbFields(columnNames.size());
  for (std::size_t i = 0; i < nbFields; i++)
  {
    vtkNew<vtkDoubleArray> arr;
    arr->SetName(columnNames[i].c_str());
    arr->SetNumberOfTuples(szX * szY);
    const std::vector<double> &data(valuesByColumn[i]);
    double *pt(arr->GetPointer(0));
    if (data.size() != szX * szY)
    {
      std::ostringstream oss;
      oss << "fillFieldInGrid : fatal internal error !"
          << "Expect " << szX << "*" << szY << "=" << szX * szY << " But having " << data.size() << " ! mail to anthony.geay@edf.fr !";
      throw MZCException(oss.str());
    }
    // transpose
    for (std::size_t i = 0; i < szX; i++)
      for (std::size_t j = 0; j < szY; j++)
        pt[j * szX + i] = data[i * szY + j];
    //
    pd->AddArray(arr);
  }
}

void vtkDepthVsTime::vtkInternal::operate(double timeStep, vtkUnstructuredGrid *usgIn)
{
  vtkNew<vtkPolyData> sourceCpy;
  sourceCpy->DeepCopy(_ladder);
  vtkNew<vtkCompositeDataProbeFilter> probeFilter;
  probeFilter->SetInputData(sourceCpy);
  probeFilter->SetSourceData(usgIn);
  probeFilter->Update();
  vtkDataObject *res(probeFilter->GetOutput());
  vtkPolyData *res2(vtkPolyData::SafeDownCast(res));
  if (!res2)
  {
    std::ostringstream oss;
    oss << "Internal error ! unexpected returned of resample filter !";
    throw MZCException(oss.str());
  }
  pushData(timeStep, res2);
}

void vtkDepthVsTime::vtkInternal::pushData(double timeStep, vtkPolyData *ds)
{
  if (!ds)
    throw MZCException("pushData : no  data !");
  vtkDataSetAttributes *dsa(ds->GetPointData());
  if (!dsa)
    throw MZCException("pushData : no point data !");
  _time.push_back(timeStep);
  if (_isInit)
    pushDataInit(timeStep, dsa);
  else
    pushDataStd(timeStep, dsa);
  _isInit = false;
}

void vtkDepthVsTime::vtkInternal::pushDataInit(double timeStep, vtkDataSetAttributes *dsa)
{
  int nba(dsa->GetNumberOfArrays());
  for (int i = 0; i < nba; i++)
  {
    vtkDataArray *arr(dsa->GetArray(i));
    if (!arr)
      continue;
    if (arr->GetNumberOfComponents() != 1)
      continue;
    std::size_t tmp(arr->GetNumberOfTuples());
    if (tmp == 0)
      continue;
    if (_nbItems == std::numeric_limits<std::size_t>::max())
      _nbItems = tmp;
    if (tmp != _nbItems)
      continue;
    const char *name(arr->GetName());
    if (!name)
      continue;
    vtkDoubleArray *arr1(vtkDoubleArray::SafeDownCast(arr));
    vtkFloatArray *arr2(vtkFloatArray::SafeDownCast(arr));
    if (!arr1 && !arr2)
      continue;
    _columnNames.push_back(name);
    if (arr1)
    {
      const double *pt(arr1->GetPointer(0));
      _data.resize(_columnNames.size());
      std::vector<double> &data(_data[_columnNames.size() - 1]);
      data.insert(data.end(), pt, pt + _nbItems);
      continue;
    }
    if (arr2)
    {
      const float *pt(arr2->GetPointer(0));
      _data.resize(_columnNames.size());
      std::vector<double> &data(_data[_columnNames.size() - 1]);
      data.insert(data.end(), pt, pt + _nbItems);
      continue;
    }
  }
}

void vtkDepthVsTime::vtkInternal::pushDataStd(double timeStep, vtkDataSetAttributes *dsa)
{
  std::set<std::string> cnsRef(_columnNames.begin(), _columnNames.end()), cns;
  int nba(dsa->GetNumberOfArrays());
  for (int i = 0; i < nba; i++)
  {
    vtkDataArray *arr(dsa->GetArray(i));
    if (!arr)
      continue;
    if (arr->GetNumberOfComponents() != 1)
      continue;
    if (arr->GetNumberOfTuples() != _nbItems)
      continue;
    const char *name(arr->GetName());
    if (!name)
      continue;
    vtkDoubleArray *arr1(vtkDoubleArray::SafeDownCast(arr));
    vtkFloatArray *arr2(vtkFloatArray::SafeDownCast(arr));
    if (!arr1 && !arr2)
      continue;
    std::string nameCpp(name);
    std::vector<std::string>::iterator it(std::find(_columnNames.begin(), _columnNames.end(), nameCpp));
    if (it == _columnNames.end())
      continue;
    std::size_t columnId(std::distance(_columnNames.begin(), it));
    if (cns.find(nameCpp) != cns.end())
      throw MZCException("pushDataStd : internal error 1 !");
    cns.insert(nameCpp);
    std::vector<double> &data(_data[columnId]);
    if (arr1)
    {
      const double *pt(arr1->GetPointer(0));
      data.insert(data.end(), pt, pt + _nbItems);
      continue;
    }
    if (arr2)
    {
      const float *pt(arr2->GetPointer(0));
      data.insert(data.end(), pt, pt + _nbItems);
      continue;
    }
  }
  if (cnsRef != cns)
    throw MZCException("Some float arrays are not present along time !");
}

void vtkDepthVsTime::vtkInternal::fillGrid(vtkRectilinearGrid *grid) const
{
  grid->SetDimensions(_time.size(), _nbItems, 1);
  {
    vtkNew<vtkDoubleArray> arrX;
    arrX->SetNumberOfTuples(_time.size());
    std::copy(_time.begin(), _time.end(), arrX->GetPointer(0));
    grid->SetXCoordinates(arrX);
  }
  {
    vtkNew<vtkDoubleArray> arrY;
    arrY->SetNumberOfTuples(_nbItems);
    double *arrYPt(arrY->GetPointer(0));
    vtkDoubleArray *data(vtkDoubleArray::SafeDownCast(_ladder->GetPoints()->GetData()));
    if (!data)
      throw MZCException("fillGrid : internal error 1 !");
    if (data->GetNumberOfTuples() != _nbItems)
      throw MZCException("fillGrid : internal error 2 !");
    const double *pt(data->GetPointer(0));
    for (std::size_t i = 0; i < _nbItems; i++)
      arrYPt[i] = pt[3 * i + 2];
    grid->SetYCoordinates(arrY);
  }
  fillFieldInGrid(grid, _data, _time.size(), _nbItems, _columnNames);
}

std::size_t vtkDepthVsTime::vtkInternal::CheckPts(vtkPointSet *usg, vtkDataArray *&arr)
{
  if (!usg)
    throw MZCException("CheckPts : expect an unstucturedgrid !");
  vtkPoints *pts(usg->GetPoints());
  if (!pts)
    throw MZCException("CheckPts : no points in grid !");
  arr = pts->GetData();
  if (!arr)
    throw MZCException("CheckPts : no data in points in grid !");
  if (arr->GetNumberOfComponents() != 3)
    throw MZCException("CheckPts : 3D expected !");
  std::size_t nbPts(arr->GetNumberOfTuples());
  if (nbPts < 1)
    throw MZCException("CheckPts : no input point !");
  vtkDoubleArray *arr1(vtkDoubleArray::SafeDownCast(arr));
  vtkFloatArray *arr2(vtkFloatArray::SafeDownCast(arr));
  if (!arr1 && !arr2)
    throw MZCException("scanCoordsOfDS : for coords expected FLOAT32 or FLOAT64 !");
  return nbPts;
}

void fillLader(vtkPolyData *source, double xpos, double ypos, double zmin, double zmax, int nbOfDiscr, const double *&ptOnInputDiscr)
{
  vtkNew<vtkDoubleArray> arr;
  arr->SetNumberOfComponents(3);
  arr->SetNumberOfTuples(nbOfDiscr);
  double *pt(arr->GetPointer(0)), delta((zmax - zmin) / ((double)(nbOfDiscr - 1)));
  for (int i = 0; i < nbOfDiscr; i++)
  {
    pt[3 * i] = xpos;
    pt[3 * i + 1] = ypos;
    pt[3 * i + 2] = ((double)i) * delta + zmin;
  }
  ptOnInputDiscr = pt;
  vtkNew<vtkPoints> pts;
  pts->SetData(arr);
  source->SetPoints(pts);
  vtkNew<vtkCellArray> verts;
  {
    vtkNew<vtkIdTypeArray> conn;
    conn->SetNumberOfComponents(1);
    conn->SetNumberOfTuples(2 * nbOfDiscr);
    vtkIdType *pt(conn->GetPointer(0));
    for (vtkIdType i = 0; i < nbOfDiscr; i++)
    {
      pt[2 * i] = 1;
      pt[2 * i + 1] = i;
    }
    verts->SetCells(nbOfDiscr, conn);
  }
  source->SetVerts(verts);
}

void vtkDepthVsTime::vtkInternal::scanCoordsOfDS(vtkUnstructuredGrid *usg, vtkPolyData *zePoint, int nbOfExpectedDiscr)
{
  vtkDataArray *arr(0);
  std::size_t nbPts(CheckPts(usg, arr));
  double Zmin(std::numeric_limits<double>::max()), Zmax(std::numeric_limits<double>::max()), LocZmin(std::numeric_limits<double>::max()), LocZmax(std::numeric_limits<double>::max());
  {
    double tmp[6];
    usg->GetBounds(tmp);
    Zmin = tmp[4];
    Zmax = tmp[5];
  }
  if (Zmin == Zmax)
    throw MZCException("scanCoordsOfDS : Zmin == Zmax ! Looks bad !");
  vtkNew<vtkUnstructuredGrid> usgCpy;
  usgCpy->DeepCopy(usg);
  vtkPointData *pd(usgCpy->GetPointData());
  if (!pd)
    throw MZCException("scanCoordsOfDS : unexpected case ! send mail to anthony.geay@edf.fr !");
  int nbArr(pd->GetNumberOfArrays());
  if (nbArr >= 1)
  {
    for (int i = nbArr - 1; i >= 0; i--)
      pd->RemoveArray(i);
  }
  {
    vtkNew<vtkDoubleArray> fakeArr;
    fakeArr->SetName("a");
    std::size_t nbPts(usgCpy->GetNumberOfPoints());
    fakeArr->SetNumberOfTuples(nbPts);
    double *pt(fakeArr->GetPointer(0));
    std::fill(pt, pt + nbPts, 1.);
    pd->AddArray(fakeArr);
  }
  //
  if (zePoint->GetNumberOfPoints() != 1)
    throw MZCException("scanCoordsOfDS : source has to have exactly one point !");
  vtkDataArray *pts(zePoint->GetPoints()->GetData());
  if (!pts)
    throw MZCException("scanCoordsOfDS : internal error ! send mail to anthony.geay@edf.fr !");
  vtkDoubleArray *pts1(vtkDoubleArray::SafeDownCast(pts));
  vtkFloatArray *pts2(vtkFloatArray::SafeDownCast(pts));
  if (!pts1 && !pts2)
    throw MZCException("scanCoordsOfDS : internal error 2 ! send mail to anthony.geay@edf.fr !");
  double Xpos(std::numeric_limits<double>::max()), Ypos(std::numeric_limits<double>::max());
  if (pts1)
  {
    Xpos = pts1->GetTypedComponent(0, 0);
    Ypos = pts1->GetTypedComponent(0, 1);
  }
  else
  {
    Xpos = (double)pts2->GetTypedComponent(0, 0);
    Ypos = (double)pts2->GetTypedComponent(0, 1);
  }
  //
  vtkNew<vtkPolyData> source;
  const double *ptOnInputDiscr(NULL);
  {
    fillLader(source, Xpos, Ypos, Zmin, Zmax, NB_OF_DISCR_TO_DEDUCE_START_STOP, ptOnInputDiscr);
  }
  //
  vtkNew<vtkCompositeDataProbeFilter> probeFilter;
  probeFilter->SetInputData(source);
  probeFilter->SetSourceData(usgCpy);
  probeFilter->Update();
  vtkDataObject *res(probeFilter->GetOutput());
  vtkPolyData *res2(vtkPolyData::SafeDownCast(res));
  if (!res2)
    throw MZCException("scanCoordsOfDS : Internal error ! unexpected returned of resample filter !");
  //
  {
    vtkPointData *pd(res2->GetPointData());
    if (!pd)
      throw MZCException("scanCoordsOfDS : internal error 3 ! send mail to anthony.geay@edf.fr !");
    vtkDataArray *pts(pd->GetArray("a"));
    if (!pts)
      throw MZCException("scanCoordsOfDS : internal error 4 ! send mail to anthony.geay@edf.fr !");
    vtkDoubleArray *pts1(vtkDoubleArray::SafeDownCast(pts));
    if (!pts1)
      throw MZCException("scanCoordsOfDS : internal error 5 ! send mail to anthony.geay@edf.fr !");
    const double *vals(pts1->GetPointer(0));
    int is(std::numeric_limits<int>::max()), ie(std::numeric_limits<int>::max());
    for (int i = 0; i < NB_OF_DISCR_TO_DEDUCE_START_STOP; i++)
    {
      if (vals[i] == 1.)
      {
        is = i;
        break;
      }
    }
    if (is == std::numeric_limits<int>::max())
      throw MZCException("scanCoordsOfDS : selected point seems to be outside of domain !");
    for (int i = NB_OF_DISCR_TO_DEDUCE_START_STOP - 1; i >= 0; i--)
    {
      if (vals[i] == 1.)
      {
        ie = i;
        break;
      }
    }
    if (is == ie)
      throw MZCException("scanCoordsOfDS : internal error 6 ! send mail to anthony.geay@edf.fr !");
    LocZmin = ptOnInputDiscr[3 * is + 2];
    LocZmax = ptOnInputDiscr[3 * ie + 2];
  }
  //
  //std::cerr << "mmmmmmm " << Xpos << " " << Ypos << std::endl;
  //std::cerr << "-> " << Zmin << " " << Zmax << std::endl;
  //std::cerr << "-> " << LocZmin << " " << LocZmax << std::endl;
  //
  {
    _ladder.TakeReference(vtkPolyData::New());
    fillLader(_ladder, Xpos, Ypos, LocZmin, LocZmax, nbOfExpectedDiscr, ptOnInputDiscr);
  }
}

int vtkDepthVsTime::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkDepthVsTime::RequestData        ##########################################" << std::endl;
  try
  {
    //
    if (this->NumberOfTimeSteps == 0)
    {
      vtkErrorMacro("No time steps in input data!");
      return 0;
    }
    vtkInformation *outInfo(outputVector->GetInformationObject(0));
    vtkUnstructuredGrid *usgIn(0);
    ExtractInfo(inputVector[0], usgIn);
    vtkInformation *sourceInfo(inputVector[1]->GetInformationObject(0));
    vtkDataObject *source(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *source2(vtkPolyData::SafeDownCast(source));
    if (!source2)
      throw MZCException("vtkPolyData expected as source !");
    // is this the first request
    if (!this->IsExecuting)
    {
      request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
      this->IsExecuting = true;
      delete this->Internal;
      this->Internal = new vtkInternal;
      this->Internal->scanCoordsOfDS(usgIn, source2, this->NbDiscrPtsAlongZ);
    }
    //
    // do something
    {
      double timeStep;
      {
        vtkInformation *inInfo(inputVector[0]->GetInformationObject(0));
        vtkDataObject *input(vtkDataObject::GetData(inInfo));
        timeStep = input->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());
      }
      this->Internal->operate(timeStep, usgIn);
    }
    this->UpdateProgress(double(this->CurrentTimeIndex) / double(this->NumberOfTimeSteps));
    //
    this->CurrentTimeIndex++;
    if (this->CurrentTimeIndex == this->NumberOfTimeSteps)
    {
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      this->CurrentTimeIndex = 0;
      this->IsExecuting = false;
      vtkInformation *outInfo(outputVector->GetInformationObject(0));
      vtkRectilinearGrid *output(vtkRectilinearGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
      vtkNew<vtkRectilinearGrid> grid;
      this->Internal->fillGrid(grid);
      output->ShallowCopy(grid);
    }
  }
  catch (MZCException &e)
  {
    if (this->IsExecuting)
    {
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      this->CurrentTimeIndex = 0;
      this->IsExecuting = false;
    }
    vtkErrorMacro(<< "Exception has been thrown in vtkDepthVsTime::RequestData : " << e.what());
    return 0;
  }
  return 1;
}

void vtkDepthVsTime::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkDepthVsTime::SetSourceData(vtkDataObject *input)
{
  this->SetInputData(1, input);
}

void vtkDepthVsTime::SetSourceConnection(vtkAlgorithmOutput *algOutput)
{
  this->SetInputConnection(1, algOutput);
}

int vtkDepthVsTime::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkRectilinearGrid");
  return 1;
}
