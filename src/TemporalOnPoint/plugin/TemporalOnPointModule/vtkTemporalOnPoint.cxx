// Copyright (C) 2021-2023  CEA, EDF
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

#include "vtkTemporalOnPoint.h"

#include <vtkAdjacentVertexIterator.h>
#include <vtkAlgorithmOutput.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCompositeDataProbeFilter.h>
#include <vtkDataArraySelection.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTable.h>
#include <vtkTimeStamp.h>
#include <vtkUnstructuredGrid.h>

#include <deque>
#include <map>
#include <set>
#include <sstream>

vtkStandardNewMacro(vtkTemporalOnPoint);

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

class vtkTemporalOnPoint::vtkInternal
{
public:
  vtkInternal()
    : _isInit(true)
  {
    _Z = std::numeric_limits<double>::max();
  }
  void operate(double timeStep, vtkUnstructuredGrid* usgIn, vtkPolyData* source);
  void pushData(double timeStep, vtkPolyData* data);
  void fillTable(vtkTable* table) const;
  void scanCoordsOfDS(vtkUnstructuredGrid* usg);
  static std::size_t CheckPts(vtkPointSet* usg, vtkDataArray*& arr);

private:
  void pushDataInit(double timeStep, vtkDataSetAttributes* dsa);
  void pushDataStd(double timeStep, vtkDataSetAttributes* dsa);

private:
  double _Z;
  bool _isInit;
  std::vector<std::string> _columnNames;
  std::vector<double> _time;
  // First level of _data is for curves series.
  // Second level of _data is for curve. Foreach i sizeof(_data[i]) must be equal to
  // sizeof(_columName) Third level of _data is for time. Foreach i,j sizeof(_data[i][j]) must be
  // equal to sizeof(_time)
  std::vector<std::vector<std::vector<double> > > _data;
};

void ExtractInfo(vtkInformationVector* inputVector, vtkUnstructuredGrid*& usgIn)
{
  vtkInformation* inputInfo(inputVector->GetInformationObject(0));
  vtkDataSet* input(0);
  vtkDataSet* input0(vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  vtkMultiBlockDataSet* input1(
    vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  if (input0)
    input = input0;
  else
  {
    if (!input1)
      throw MZCException(
        "Input dataSet must be a DataSet or single elt multi block dataset expected !");
    if (input1->GetNumberOfBlocks() != 1)
      throw MZCException("Input dataSet is a multiblock dataset with not exactly one block ! Use "
                         "MergeBlocks or ExtractBlocks filter before calling this filter !");
    vtkDataObject* input2(input1->GetBlock(0));
    if (!input2)
      throw MZCException("Input dataSet is a multiblock dataset with exactly one block but this "
                         "single element is NULL !");
    vtkDataSet* input2c(vtkDataSet::SafeDownCast(input2));
    if (!input2c)
      throw MZCException(
        "Input dataSet is a multiblock dataset with exactly one block but this single element is "
        "not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    input = input2c;
  }
  if (!input)
    throw MZCException("Input data set is NULL !");
  usgIn = vtkUnstructuredGrid::SafeDownCast(input);
  if (!usgIn)
    throw MZCException("Input data set is not an unstructured mesh ! This filter works only on "
                       "unstructured meshes !");
}

////////////////////

vtkTemporalOnPoint::vtkTemporalOnPoint()
  : NumberOfTimeSteps(0)
  , IsExecuting(false)
  , CurrentTimeIndex(0)
  , Internal(nullptr)
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

vtkTemporalOnPoint::~vtkTemporalOnPoint()
{
  delete this->Internal;
  this->Internal = nullptr;
}

int vtkTemporalOnPoint::RequestInformation(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // std::cerr << "########################################## vtkTemporalOnPoint::RequestInformation
  // ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid* usgIn(0);
    ExtractInfo(inputVector[0], usgIn);
    vtkInformation* inInfo(inputVector[0]->GetInformationObject(0));
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
    vtkInformation* outInfo(outputVector->GetInformationObject(0));
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
  catch (MZCException& e)
  {
    vtkErrorMacro(<< "Exception has been thrown in vtkTemporalOnPoint::RequestInformation : "
                  << e.what());
    return 0;
  }
  return 1;
}

int vtkTemporalOnPoint::RequestUpdateExtent(vtkInformation*, vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed(outputVector))
{
  // vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkInformation* inInfo1 = inputVector[0]->GetInformationObject(0);

  // get the requested update extent
  double* inTimes = inInfo1->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  if (inTimes)
  {
    double timeReq = inTimes[this->CurrentTimeIndex];
    inInfo1->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), timeReq);
  }

  return 1;
}

std::string buildNameOfEntryFrom(const std::string& name, std::size_t id, std::size_t nbOfElt)
{
  if (nbOfElt == 0)
    throw MZCException("buildNameOfEntryFrom : nbElt == 0 !");
  if (nbOfElt == 1)
    return name;
  std::ostringstream oss;
  oss << name << "_" << id;
  return oss.str();
}

void buildTableFrom(vtkTable* table,
  const std::vector<std::vector<std::vector<double> > >& valuesByColumn,
  const std::vector<double> zeTimes, const std::vector<std::string>& columnNames)
{
  {
    vtkNew<vtkDoubleArray> timeArr;
    timeArr->SetName("Time");
    timeArr->SetNumberOfTuples(zeTimes.size());
    double* pt(timeArr->GetPointer(0));
    std::copy(zeTimes.begin(), zeTimes.end(), pt);
    table->AddColumn(timeArr);
  }
  std::size_t nbOfSeries(valuesByColumn.size()), nbCols(columnNames.size());
  for (std::size_t i = 0; i < nbOfSeries; i++)
  {
    for (std::size_t j = 0; j < nbCols; j++)
    {
      vtkNew<vtkDoubleArray> arr;
      std::string name(buildNameOfEntryFrom(columnNames[j], i, nbOfSeries));
      arr->SetName(name.c_str());
      arr->SetNumberOfTuples(zeTimes.size());
      double* pt(arr->GetPointer(0));
      std::copy(valuesByColumn[i][j].begin(), valuesByColumn[i][j].end(), pt);
      table->AddColumn(arr);
    }
  }
}

void vtkTemporalOnPoint::vtkInternal::operate(
  double timeStep, vtkUnstructuredGrid* usgIn, vtkPolyData* source)
{
  vtkNew<vtkPolyData> sourceCpy;
  sourceCpy->DeepCopy(source);
  vtkDataArray* arr(0);
  std::size_t nbPts(CheckPts(sourceCpy, arr));
  if (_Z != std::numeric_limits<double>::max())
  {
    vtkDoubleArray* arr1(vtkDoubleArray::SafeDownCast(arr));
    vtkFloatArray* arr2(vtkFloatArray::SafeDownCast(arr));
    if (arr1)
    {
      double* pt(arr1->GetPointer(0));
      for (std::size_t i = 0; i < nbPts; i++)
        pt[3 * i + 2] = _Z;
    }
    else
    {
      float* pt(arr2->GetPointer(0));
      for (std::size_t i = 0; i < nbPts; i++)
        pt[3 * i + 2] = _Z;
    }
  }
  //
  vtkNew<vtkCompositeDataProbeFilter> probeFilter;
  probeFilter->SetInputData(sourceCpy);
  probeFilter->SetSourceData(usgIn);
  probeFilter->Update();
  vtkDataObject* res(probeFilter->GetOutput());
  vtkPolyData* res2(vtkPolyData::SafeDownCast(res));
  if (!res2)
  {
    std::ostringstream oss;
    oss << "Internal error ! unexpected returned of resample filter !";
    throw MZCException(oss.str());
  }
  pushData(timeStep, res2);
}

void vtkTemporalOnPoint::vtkInternal::pushData(double timeStep, vtkPolyData* ds)
{
  if (!ds)
    throw MZCException("pushData : no  data !");
  vtkDataSetAttributes* dsa(ds->GetPointData());
  if (!dsa)
    throw MZCException("pushData : no point data !");
  _time.push_back(timeStep);
  if (_isInit)
    pushDataInit(timeStep, dsa);
  else
    pushDataStd(timeStep, dsa);
  _isInit = false;
}

void vtkTemporalOnPoint::vtkInternal::pushDataInit(double timeStep, vtkDataSetAttributes* dsa)
{
  std::size_t nbOfItems(std::numeric_limits<std::size_t>::max());
  int nba(dsa->GetNumberOfArrays());
  for (int i = 0; i < nba; i++)
  {
    vtkDataArray* arr(dsa->GetArray(i));
    if (!arr)
      continue;
    if (arr->GetNumberOfComponents() != 1)
      continue;
    std::size_t tmp(arr->GetNumberOfTuples());
    if (tmp == 0)
      continue;
    if (nbOfItems == std::numeric_limits<std::size_t>::max())
    {
      nbOfItems = tmp;
      _data.resize(nbOfItems);
    }
    if (tmp != nbOfItems)
      continue;
    const char* name(arr->GetName());
    if (!name)
      continue;
    vtkDoubleArray* arr1(vtkDoubleArray::SafeDownCast(arr));
    vtkFloatArray* arr2(vtkFloatArray::SafeDownCast(arr));
    if (!arr1 && !arr2)
      continue;
    _columnNames.push_back(name);
    if (arr1)
    {
      const double* pt(arr1->GetPointer(0));
      for (std::size_t j = 0; j < nbOfItems; j++)
      {
        _data[j].resize(_columnNames.size());
        _data[j][_columnNames.size() - 1].push_back(pt[j]);
      }
      continue;
    }
    if (arr2)
    {
      const float* pt(arr2->GetPointer(0));
      for (std::size_t j = 0; j < nbOfItems; j++)
      {
        _data[j].resize(_columnNames.size());
        _data[j][_columnNames.size() - 1].push_back(pt[j]);
      }
      continue;
    }
  }
}

void vtkTemporalOnPoint::vtkInternal::pushDataStd(double timeStep, vtkDataSetAttributes* dsa)
{
  std::set<std::string> cnsRef(_columnNames.begin(), _columnNames.end()), cns;
  std::size_t nbOfItems(_data.size());
  int nba(dsa->GetNumberOfArrays());
  for (int i = 0; i < nba; i++)
  {
    vtkDataArray* arr(dsa->GetArray(i));
    if (!arr)
      continue;
    if (arr->GetNumberOfComponents() != 1)
      continue;
    if (arr->GetNumberOfTuples() != nbOfItems)
      continue;
    const char* name(arr->GetName());
    if (!name)
      continue;
    vtkDoubleArray* arr1(vtkDoubleArray::SafeDownCast(arr));
    vtkFloatArray* arr2(vtkFloatArray::SafeDownCast(arr));
    if (!arr1 && !arr2)
      continue;
    std::string nameCpp(name);
    std::vector<std::string>::iterator it(
      std::find(_columnNames.begin(), _columnNames.end(), nameCpp));
    if (it == _columnNames.end())
      continue;
    std::size_t columnId(std::distance(_columnNames.begin(), it));
    cns.insert(nameCpp);
    if (arr1)
    {
      const double* pt(arr1->GetPointer(0));
      for (std::size_t j = 0; j < nbOfItems; j++)
        _data[j][columnId].push_back(pt[j]);
      continue;
    }
    if (arr2)
    {
      const float* pt(arr2->GetPointer(0));
      for (std::size_t j = 0; j < nbOfItems; j++)
        _data[j][columnId].push_back(pt[j]);
      continue;
    }
  }
  if (cnsRef != cns)
    throw MZCException("Some float arrays are not present along time !");
}

void vtkTemporalOnPoint::vtkInternal::fillTable(vtkTable* table) const
{
  buildTableFrom(table, _data, _time, _columnNames);
}

std::size_t vtkTemporalOnPoint::vtkInternal::CheckPts(vtkPointSet* usg, vtkDataArray*& arr)
{
  if (!usg)
    throw MZCException("CheckPts : expect an unstucturedgrid !");
  vtkPoints* pts(usg->GetPoints());
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
  vtkDoubleArray* arr1(vtkDoubleArray::SafeDownCast(arr));
  vtkFloatArray* arr2(vtkFloatArray::SafeDownCast(arr));
  if (!arr1 && !arr2)
    throw MZCException("scanCoordsOfDS : for coords expected FLOAT32 or FLOAT64 !");
  return nbPts;
}

void vtkTemporalOnPoint::vtkInternal::scanCoordsOfDS(vtkUnstructuredGrid* usg)
{
  vtkDataArray* arr(0);
  std::size_t nbPts(CheckPts(usg, arr));
  vtkDoubleArray* arr1(vtkDoubleArray::SafeDownCast(arr));
  vtkFloatArray* arr2(vtkFloatArray::SafeDownCast(arr));
  if (arr1)
  {
    const double* pt(arr1->GetPointer(0));
    double ref(pt[2]);
    for (std::size_t i = 1; i < nbPts; i++)
      if (pt[3 * i + 2] != ref)
      {
        _Z = std::numeric_limits<double>::max();
        return;
      }
    _Z = ref;
  }
  else
  {
    const float* pt(arr2->GetPointer(0));
    float ref(pt[2]);
    for (std::size_t i = 1; i < nbPts; i++)
      if (pt[3 * i + 2] != ref)
      {
        _Z = std::numeric_limits<double>::max();
        return;
      }
    _Z = ref;
  }
  // std::cerr << "Is 2D ?  " << _Z << std::endl;
}

int vtkTemporalOnPoint::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // std::cerr << "########################################## vtkTemporalOnPoint::RequestData
  // ##########################################" << std::endl;
  try
  {
    //
    if (this->NumberOfTimeSteps == 0)
    {
      vtkErrorMacro("No time steps in input data!");
      return 0;
    }
    vtkInformation* outInfo(outputVector->GetInformationObject(0));
    vtkUnstructuredGrid* usgIn(0);
    ExtractInfo(inputVector[0], usgIn);
    // is this the first request
    if (!this->IsExecuting)
    {
      request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
      this->IsExecuting = true;
      delete this->Internal;
      this->Internal = new vtkInternal;
      this->Internal->scanCoordsOfDS(usgIn);
    }
    //
    // do something
    {
      vtkInformation* sourceInfo(inputVector[1]->GetInformationObject(0));
      vtkDataObject* source(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
      vtkPolyData* source2(vtkPolyData::SafeDownCast(source));
      if (!source2)
        throw MZCException("vtkPolyData expected as source !");
      double timeStep;
      {
        vtkInformation* inInfo(inputVector[0]->GetInformationObject(0));
        vtkDataObject* input(vtkDataObject::GetData(inInfo));
        timeStep = input->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());
      }
      this->Internal->operate(timeStep, usgIn, source2);
    }
    //
    this->CurrentTimeIndex++;
    if (this->CurrentTimeIndex == this->NumberOfTimeSteps)
    {
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      this->CurrentTimeIndex = 0;
      this->IsExecuting = false;
      vtkInformation* outInfo(outputVector->GetInformationObject(0));
      vtkTable* output(vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
      vtkNew<vtkTable> table;
      this->Internal->fillTable(table);
      output->ShallowCopy(table);
    }
  }
  catch (MZCException& e)
  {
    if (this->IsExecuting)
    {
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      this->CurrentTimeIndex = 0;
      this->IsExecuting = false;
    }
    vtkErrorMacro(<< "Exception has been thrown in vtkTemporalOnPoint::RequestData : " << e.what());
    return 0;
  }
  return 1;
}

void vtkTemporalOnPoint::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkTemporalOnPoint::SetSourceData(vtkDataObject* input)
{
  this->SetInputData(1, input);
}

void vtkTemporalOnPoint::SetSourceConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}

int vtkTemporalOnPoint::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
  return 1;
}
