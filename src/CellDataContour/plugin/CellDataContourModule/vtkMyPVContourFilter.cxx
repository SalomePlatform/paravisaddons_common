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

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMyContourFilter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkMyPVContourFilter.h"

#include <vtkAMRDualContour.h>
#include <vtkAppendPolyData.h>
#include <vtkCompositeDataIterator.h>
#include <vtkDataObject.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkHierarchicalBoxDataSet.h>
#include <vtkInformation.h>
#include <vtkInformationStringVectorKey.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>

vtkStandardNewMacro(vtkMyPVContourFilter);

//-----------------------------------------------------------------------------
vtkMyPVContourFilter::vtkMyPVContourFilter()
  : vtkMyContourFilter()
{
  this->SetComputeNormals(true);
}

//-----------------------------------------------------------------------------
vtkMyPVContourFilter::~vtkMyPVContourFilter()
{
}

//-----------------------------------------------------------------------------
void vtkMyPVContourFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//-----------------------------------------------------------------------------
int vtkMyPVContourFilter::ProcessRequest(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // create the output
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
  {
    return this->RequestDataObject(request, inputVector, outputVector);
  }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//-----------------------------------------------------------------------------
int vtkMyPVContourFilter::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  if (!inInfo)
  {
    vtkErrorMacro("Failed to get input information.");
    return 1;
  }

  vtkDataObject* inDataObj = inInfo->Get(vtkDataObject::DATA_OBJECT());
  if (!inDataObj)
  {
    vtkErrorMacro("Failed to get input data object.");
    return 1;
  }

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  if (!outInfo)
  {
    vtkErrorMacro("Failed to get output information.");
    return 1;
  }

  vtkDataObject* outDataObj = outInfo->Get(vtkDataObject::DATA_OBJECT());
  if (!outDataObj)
  {
    vtkErrorMacro("Failed get output data object.");
    return 1;
  }

  // Check if input is AMR data.
  if (vtkHierarchicalBoxDataSet::SafeDownCast(inDataObj))
  {
    // This is a lot to go through to get the name of the array to process.
    vtkInformation* inArrayInfo = this->GetInputArrayInformation(0);
    if (!inArrayInfo)
    {
      vtkErrorMacro("Problem getting name of array to process.");
      return 0;
    }
    int fieldAssociation = -1;
    if (!inArrayInfo->Has(vtkDataObject::FIELD_ASSOCIATION()))
    {
      vtkErrorMacro("Unable to query field association for the scalar.");
      return 0;
    }
    fieldAssociation = inArrayInfo->Get(vtkDataObject::FIELD_ASSOCIATION());
    if (fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_CELLS)
    {
      vtkSmartPointer<vtkAMRDualContour> amrDC(vtkSmartPointer<vtkAMRDualContour>::New());

      amrDC->SetInputData(0, inDataObj);
      amrDC->SetInputArrayToProcess(0, inArrayInfo);
      amrDC->SetEnableCapping(1);
      amrDC->SetEnableDegenerateCells(1);
      amrDC->SetEnableMultiProcessCommunication(1);
      amrDC->SetSkipGhostCopy(1);
      amrDC->SetTriangulateCap(1);
      amrDC->SetEnableMergePoints(1);

      for (int i = 0; i < this->GetNumberOfContours(); ++i)
      {
        vtkSmartPointer<vtkMultiBlockDataSet> out(vtkSmartPointer<vtkMultiBlockDataSet>::New());
        amrDC->SetIsoValue(this->GetValue(i));
        amrDC->Update();
        out->ShallowCopy(amrDC->GetOutput(0));
        vtkMultiBlockDataSet::SafeDownCast(outDataObj)->SetBlock(i, out);
      }
      return 1;
    }
  }

  return this->ContourUsingSuperclass(request, inputVector, outputVector);
}

//-----------------------------------------------------------------------------
int vtkMyPVContourFilter::RequestDataObject(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  if (!inInfo)
  {
    return 0;
  }

  vtkHierarchicalBoxDataSet* input = vtkHierarchicalBoxDataSet::GetData(inInfo);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  if (input)
  {
    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::GetData(outInfo);
    if (!output)
    {
      output = vtkMultiBlockDataSet::New();
      outInfo->Set(vtkDataObject::DATA_OBJECT(), output);
      this->GetOutputPortInformation(0)->Set(
        vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType());
      output->Delete();
    }
    return 1;
  }
  else
  {
    vtkDataSet* output = vtkDataSet::GetData(outInfo);
    if (!output)
    {
      output = vtkPolyData::New();
      outInfo->Set(vtkDataObject::DATA_OBJECT(), output);
      this->GetOutputPortInformation(0)->Set(
        vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType());
      output->Delete();
    }
    return 1;
  }
}

//----------------------------------------------------------------------------
int vtkMyPVContourFilter::ContourUsingSuperclass(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkDataObject* inputDO = vtkDataObject::GetData(inputVector[0], 0);
  vtkDataObject* outputDO = vtkDataObject::GetData(outputVector, 0);

  vtkCompositeDataSet* inputCD = vtkCompositeDataSet::SafeDownCast(inputDO);
  if (!inputCD)
  {
    return this->Superclass::RequestData(request, inputVector, outputVector);
  }

  vtkCompositeDataSet* outputCD = vtkCompositeDataSet::SafeDownCast(outputDO);
  outputCD->CopyStructure(inputCD);

  vtkSmartPointer<vtkCompositeDataIterator> iter;
  iter.TakeReference(inputCD->NewIterator());

  // for input.
  vtkSmartPointer<vtkInformationVector> newInInfoVec = vtkSmartPointer<vtkInformationVector>::New();
  vtkSmartPointer<vtkInformation> newInInfo = vtkSmartPointer<vtkInformation>::New();
  newInInfoVec->SetInformationObject(0, newInInfo);

  // for output.
  vtkSmartPointer<vtkInformationVector> newOutInfoVec =
    vtkSmartPointer<vtkInformationVector>::New();
  vtkSmartPointer<vtkInformation> newOutInfo = vtkSmartPointer<vtkInformation>::New();
  newOutInfoVec->SetInformationObject(0, newOutInfo);

  // Loop over all the datasets.
  for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
  {
    newInInfo->Set(vtkDataObject::DATA_OBJECT(), iter->GetCurrentDataObject());
    vtkPolyData* polydata = vtkPolyData::New();
    newOutInfo->Set(vtkDataObject::DATA_OBJECT(), polydata);
    polydata->FastDelete();

    vtkInformationVector* newInInfoVecPtr = newInInfoVec.GetPointer();
    if (!this->Superclass::RequestData(request, &newInInfoVecPtr, newOutInfoVec.GetPointer()))
    {
      return 0;
    }
    outputCD->SetDataSet(iter, polydata);
  }

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMyPVContourFilter::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
  return 1;
}

//-----------------------------------------------------------------------------
int vtkMyPVContourFilter::FillInputPortInformation(int port, vtkInformation* info)
{
  this->Superclass::FillInputPortInformation(port, info);

  // According to the documentation this is the way to append additional
  // input data set type since VTK 5.2.
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkHierarchicalBoxDataSet");
  return 1;
}