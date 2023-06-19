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

/*=========================================================================

  Program:   ParaView
  Module:    vtkAppendAttributesOverTime.cxx

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkAppendAttributesOverTime.h"

#include <vtkCompositeDataIterator.h>
#include <vtkCompositeDataSet.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkFieldData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <cmath>
#include <sstream>
#include <iomanip>

vtkStandardNewMacro(vtkAppendAttributesOverTime);

//----------------------------------------------------------------------------
vtkAppendAttributesOverTime::vtkAppendAttributesOverTime()
  :RestoreOriginalTimeStep(false)
{
}

//----------------------------------------------------------------------------
namespace
{
std::string GetPaddedValue(int value, int maxValue)
{
  int padding = 1 + static_cast<int>(log10(maxValue));
  std::ostringstream os;
  os << std::setw(padding) << std::setfill('0') << std::to_string(value);
  return os.str();
}

//----------------------------------------------------------------------------
void ResetAttributes(vtkDataObject* obj)
{
  vtkCompositeDataSet* cObj = vtkCompositeDataSet::SafeDownCast(obj);
  if (cObj)
  {
    vtkSmartPointer<vtkCompositeDataIterator> iter;
    iter.TakeReference(cObj->NewIterator());
    iter->InitTraversal();
    for (; !iter->IsDoneWithTraversal(); iter->GoToNextItem())
    {
      ResetAttributes(vtkDataSet::SafeDownCast(iter->GetCurrentDataObject()));
    }
  }
  else
  {
    // reset fields data, so TempDataObject contains only geometry.
    for (int attr = 0; attr < vtkDataObject::NUMBER_OF_ATTRIBUTE_TYPES; attr++)
    {
      vtkFieldData *fd = obj->GetAttributesAsFieldData(attr);
      if (fd != nullptr)
      {
        fd->Initialize();
      }
    }
  }
}
}

//----------------------------------------------------------------------------
bool vtkAppendAttributesOverTime::GetOutputArrayName(vtkFieldData* vtkNotUsed(arrays),
  const char* arrayName, int vtkNotUsed(inputIndex), std::string& outArrayName)
{
  std::string inString = ::GetPaddedValue(this->CurrentInputIndex, this->TimeSteps.size());
  std::string tsString = ::GetPaddedValue(this->UpdateTimeIndex, this->TimeSteps[this->CurrentInputIndex].size());
  outArrayName = std::string(arrayName) + "_input_" + inString + "_ts_" + tsString;
  return true;
}

//----------------------------------------------------------------------------
int vtkAppendAttributesOverTime::RequestUpdateExtent(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inInfo, vtkInformationVector* vtkNotUsed(outInfo))
{
  vtkInformation* info = inInfo[0]->GetInformationObject(this->CurrentInputIndex);
  if (this->RestoreOriginalTimeStep)
  {
    info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->OriginalTimeStep);
  }
  else
  {
    this->OriginalTimeStep = info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(),
              this->TimeSteps[this->CurrentInputIndex][this->UpdateTimeIndex]);
  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkAppendAttributesOverTime::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inInfo, vtkInformationVector* outInfoVec)
{
  this->UpdateTimeIndex = 0;
  this->CurrentInputIndex = 0;
  int num = inInfo[0]->GetNumberOfInformationObjects();

  // Fill this->TimeSteps
  for (int idx = 0; idx < num; idx++)
  {
    vtkInformation* info = inInfo[0]->GetInformationObject(idx);
    int len = info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    double* timeSteps = info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    std::vector<double> timeStepsVector;
    timeStepsVector.resize(len);
    std::copy(timeSteps, timeSteps + len, timeStepsVector.begin());
    this->TimeSteps.push_back(timeStepsVector);
  }

  vtkInformation* outInfo = outInfoVec->GetInformationObject(0);
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

  return 1;
}

//----------------------------------------------------------------------------
/**
 * This method calls the Superclass::RequestData to perform the real
 * merge.
 * Here we do the following:
 * - extract a dataset from 2 internal vars: an input index and a timestep index.
 * (see this->TimeSteps).
 * - call the Superclass RequestData method with 2 datasets as inputs:
 * buffer containing previously merged arrays and current dataset to process.
 * - store the resulting dataset in the buffer
 * - increment the timestep index. When last timestep of current input was
 * processed, increment input index and reset timestep index
 * - loop
 */
int vtkAppendAttributesOverTime::RequestData(
  vtkInformation* request, vtkInformationVector** inInfo, vtkInformationVector* outInfoVec)
{
  if (this->RestoreOriginalTimeStep)
  {
    this->RestoreOriginalTimeStep = false;
    request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
    return 1;
  }

  vtkInformation* info = inInfo[0]->GetInformationObject(this->CurrentInputIndex);
  vtkInformation* outInfo = outInfoVec->GetInformationObject(0);

  if (this->UpdateTimeIndex == 0 && this->CurrentInputIndex == 0)
  {
    // Before looping, we create the temporary output object
    vtkInformation *inputInfo = inInfo[0]->GetInformationObject(0);
    vtkDataObject *inputData = vtkDataObject::GetData(inputInfo);
    this->TempDataObject = vtkSmartPointer<vtkDataObject>::Take(inputData->NewInstance());
    this->TempDataObject->ShallowCopy(inputData);
    // Remove all arrays in the ouput
    ::ResetAttributes(this->TempDataObject);
  }

  this->CurrentOutInfo->Set(vtkDataObject::DATA_OBJECT(), this->TempDataObject);
  vtkInformationVector* reducedInputVec = vtkInformationVector::New();
  reducedInputVec->Append(this->CurrentOutInfo);
  reducedInputVec->Append(info);

  // perform the effective merge.
  this->Superclass::RequestData(request, &reducedInputVec, outInfoVec);
  reducedInputVec->Delete();
  // save merge result in a buffer.
  this->TempDataObject->ShallowCopy(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (this->UpdateTimeIndex <
    static_cast<vtkIdType>(this->TimeSteps.at(this->CurrentInputIndex).size()) - 1)
  {
    this->UpdateTimeIndex++;
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
  }
  else if (this->CurrentInputIndex < static_cast<vtkIdType>(this->TimeSteps.size()) - 1)
  {
    this->UpdateTimeIndex = 0;
    this->CurrentInputIndex++;
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
  }
  else
  {
    this->RestoreOriginalTimeStep = true;
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
    vtkDataObject* output = outInfo->Get(vtkDataObject::DATA_OBJECT());
    output->ShallowCopy(this->TempDataObject);
  }

  return 1;
}
