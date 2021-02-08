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

#include "vtkMoveZCote.h"

#include <vtkAdjacentVertexIterator.h>
#include <vtkAlgorithmOutput.h>
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArraySelection.h>
#include <vtkDataObjectTreeIterator.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkInEdgeIterator.h>
#include <vtkInformation.h>
#include <vtkInformationDataObjectKey.h>
#include <vtkInformationStringKey.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkTimeStamp.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVariantArray.h>
#include <vtkWarpScalar.h>

#include <deque>
#include <map>
#include <sstream>

vtkStandardNewMacro(vtkMoveZCote);

static const char ZE_ARRAY_NAME[] = "COTE Z";

static const char ZE_DISPLACEMENT_NAME[] = "@@DisplacementZ?@@";

///////////////////

class MZCException : public std::exception
{
public:
  MZCException(const std::string &s) : _reason(s) {}
  virtual const char *what() const throw() { return _reason.c_str(); }
  virtual ~MZCException() throw() {}

private:
  std::string _reason;
};

void ExtractInfo(vtkInformationVector *inputVector, vtkUnstructuredGrid *&usgIn, vtkDoubleArray *&arr)
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
  vtkPointData *att(usgIn->GetPointData());
  if (!att)
    throw MZCException("Input dataset has no point data attribute ! Impossible to move mesh !");
  vtkDataArray *zeArr(0);
  for (int i = 0; i < att->GetNumberOfArrays(); i++)
  {
    vtkDataArray *locArr(att->GetArray(i));
    std::string s(locArr->GetName());
    if (s == ZE_ARRAY_NAME)
    {
      zeArr = locArr;
      break;
    }
  }
  if (!zeArr)
  {
    std::ostringstream oss;
    oss << "Impossible to locate the array called \"" << ZE_ARRAY_NAME << "\" used to move mesh !";
    throw MZCException(oss.str());
  }
  arr = vtkDoubleArray::SafeDownCast(zeArr);
  if (!arr)
  {
    std::ostringstream oss;
    oss << "Array called \"" << ZE_ARRAY_NAME << "\" has been located but this is NOT a float64 array !";
    throw MZCException(oss.str());
  }
  if (arr->GetNumberOfComponents() != 1)
  {
    std::ostringstream oss;
    oss << "Float64 array called \"" << ZE_ARRAY_NAME << "\" has been located but this array has not exactly 1 components as it should !";
    throw MZCException(oss.str());
  }
  if (arr->GetNumberOfTuples() != usgIn->GetNumberOfPoints())
  {
    std::ostringstream oss;
    oss << "Float64-1 components array called \"" << ZE_ARRAY_NAME << "\" has been located but the number of tuples is invalid ! Should be " << usgIn->GetNumberOfPoints() << " instead of " << arr->GetNumberOfTuples() << " !";
    throw MZCException(oss.str());
  }
}

////////////////////
int vtkMoveZCote::RequestInformation(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkMoveZCote::RequestInformation ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid *usgIn(0);
    vtkDoubleArray *arr(0);
    ExtractInfo(inputVector[0], usgIn, arr);
  }
  catch (MZCException &e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkMoveZCote::RequestInformation : " << e.what() << std::endl;
    if (this->HasObserver("ErrorEvent"))
    {
      this->InvokeEvent("ErrorEvent", const_cast<char *>(oss.str().c_str()));
    }
    else
    {
      vtkOutputWindowDisplayErrorText(const_cast<char *>(oss.str().c_str()));
    }
    vtkObject::BreakOnError();
    return 0;
  }
  return 1;
}

int vtkMoveZCote::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkMoveZCote::RequestData        ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid *usgIn(0);
    vtkDoubleArray *arr(0);
    ExtractInfo(inputVector[0], usgIn, arr);
    //
    int nbPts(usgIn->GetNumberOfPoints());
    vtkSmartPointer<vtkUnstructuredGrid> step1(vtkSmartPointer<vtkUnstructuredGrid>::New());
    step1->DeepCopy(usgIn);
    vtkSmartPointer<vtkDoubleArray> arr1(vtkSmartPointer<vtkDoubleArray>::New());
    arr1->SetName(ZE_DISPLACEMENT_NAME);
    arr1->SetNumberOfComponents(1);
    arr1->SetNumberOfTuples(nbPts);
    double *ptToFeed(arr1->Begin());
    vtkDataArray *coords(usgIn->GetPoints()->GetData());
    vtkDoubleArray *coords2(vtkDoubleArray::SafeDownCast(coords));
    if (!coords2)
    {
      throw MZCException("Input coordinates are not float64 !");
    }
    if (coords2->GetNumberOfComponents() != 3)
    {
      throw MZCException("Input coordinates do not have 3 components as it should !");
    }
    const double *srcPt1(arr->Begin()), *srcPt2(coords2->Begin());
    for (int i = 0; i < nbPts; i++, ptToFeed++)
    {
      *ptToFeed = srcPt1[i] - srcPt2[3 * i + 2];
    }
    int idx(step1->GetPointData()->AddArray(arr1));
    step1->GetPointData()->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    //
    vtkSmartPointer<vtkWarpScalar> ws(vtkSmartPointer<vtkWarpScalar>::New());
    ws->SetInputData(step1);
    ws->SetScaleFactor(1);
    ws->SetInputArrayToProcess(idx, 0, 0, "vtkDataObject::FIELD_ASSOCIATION_POINTS", ZE_DISPLACEMENT_NAME);
    ws->Update();
    vtkSmartPointer<vtkDataSet> ds(ws->GetOutput());
    //
    ds->GetPointData()->RemoveArray(idx);
    //
    vtkInformation *outInfo(outputVector->GetInformationObject(0));
    vtkUnstructuredGrid *output(vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
    //
    output->ShallowCopy(ds);
  }
  catch (MZCException &e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkMoveZCote::RequestInformation : " << e.what() << std::endl;
    if (this->HasObserver("ErrorEvent"))
    {
      this->InvokeEvent("ErrorEvent", const_cast<char *>(oss.str().c_str()));
    }
    else
    {
      vtkOutputWindowDisplayErrorText(const_cast<char *>(oss.str().c_str()));
    }
    vtkObject::BreakOnError();
    return 0;
  }
  return 1;
}

void vtkMoveZCote::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
