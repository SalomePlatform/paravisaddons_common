// Copyright (C) 2021-2022  CEA/DEN, EDF R&D
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

#include "vtkElectromagnetismStreamTraceur.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkMergeBlocks.h>
#include <vtkMultiBlockDataGroupFilter.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>
#include "vtkStreamTracer.h"
#include "vtkPointSource.h"
#include "vtkPCellDataToPointData.h"

#include "vtkCompositeDataIterator.h"
#include "vtkCompositeInterpolatedVelocityField.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkInterpolatedVelocityField.h"

#include <vector>

vtkObjectFactoryNewMacro(vtkElectromagnetismStreamTraceur);

vtkElectromagnetismStreamTraceur::vtkElectromagnetismStreamTraceur():IntegrationDirection(BOTH),IntegratorType(RUNGE_KUTTA45),IntegrationStepUnit(2)
,InitialIntegrationStep(0.2),TerminalSpeed(1e-12),MaximumError(1e-6),MaximumNumberOfSteps(2000)
{
  this->SetNumberOfInputPorts(2);
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::VECTORS);
}

void vtkElectromagnetismStreamTraceur::SetSourceConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}

void vtkElectromagnetismStreamTraceur::SetSourceData(vtkDataSet* source)
{
  this->SetInputData(1, source);
}

vtkDataSet* vtkElectromagnetismStreamTraceur::GetSource()
{
  if (this->GetNumberOfInputConnections(1) < 1)
  {
    return nullptr;
  }
  return vtkDataSet::SafeDownCast(this->GetExecutive()->GetInputData(1, 0));
}

int vtkElectromagnetismStreamTraceur::SetupOutput(vtkInformation* inInfo, vtkInformation* outInfo)
{
  int piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  vtkDataObject* input = inInfo->Get(vtkDataObject::DATA_OBJECT());
  vtkDataObject* output = outInfo->Get(vtkDataObject::DATA_OBJECT());

  // Pass through field data
  output->GetFieldData()->PassData(input->GetFieldData());

  vtkCompositeDataSet* hdInput = vtkCompositeDataSet::SafeDownCast(input);
  vtkDataSet* dsInput = vtkDataSet::SafeDownCast(input);
  if (hdInput)
  {
    this->InputData = hdInput;
    hdInput->Register(this);
    return 1;
  }
  else if (dsInput)
  {
    vtkNew<vtkMultiBlockDataSet> mb;
    mb->SetNumberOfBlocks(numPieces);
    mb->SetBlock(piece, dsInput);
    this->InputData = mb;
    mb->Register(this);
    return 1;
  }
  else
  {
    vtkErrorMacro(
      "This filter cannot handle input of type: " << (input ? input->GetClassName() : "(none)"));
    return 0;
  }
}

int vtkElectromagnetismStreamTraceur::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkDataArray* arr = nullptr;
  vtkDataSet* input0 = nullptr;
  if (!this->SetupOutput(inInfo, outInfo))
  {
    return 0;
  }

  vtkInformation* sourceInfo = inputVector[1]->GetInformationObject(0);
  vtkDataSet* source = nullptr;
  if (sourceInfo)
  {
    source = vtkDataSet::SafeDownCast(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
  }

  vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkSmartPointer<vtkCompositeDataIterator> iterP;
  iterP.TakeReference(this->InputData->NewIterator());

  iterP->GoToFirstItem();
  if (!iterP->IsDoneWithTraversal() && !input0)
  {
    input0 = vtkDataSet::SafeDownCast(iterP->GetCurrentDataObject());
    iterP->GoToNextItem();
  }

  int vecType(0);
  arr =  this->GetInputArrayToProcess(0, input0, vecType);
  if(!arr)
  {
    vtkErrorMacro("No vector field selected in input !");
    return 0;
  }

  vtkDataSet *input( vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT())) );
  const char *ArrayForGlyph(arr->GetName());
  //
  vtkNew<vtkPCellDataToPointData> cc;
  cc->SetInputData(input0);
  cc->SetProcessAllArrays(1);
  cc->SetPassCellData(0);
  cc->SetPieceInvariant(0);
  cc->Update();
  // Compute default Maximum Propagation
  input0->ComputeBounds();
  double james[6];
  input0->GetBounds(james);
  double dftMaxProp(std::min(std::min(james[1]-james[0],james[3]-james[2]),james[5]-james[4]));
  //
  vtkNew<vtkStreamTracer> streamTracer;
  streamTracer->SetInputConnection(cc->GetOutputPort());
  streamTracer->SetInterpolatorTypeToDataSetPointLocator();
  streamTracer->SetIntegrationDirection(this->IntegrationDirection);
  streamTracer->SetIntegratorType(this->IntegratorType);
  streamTracer->SetIntegrationStepUnit(this->IntegrationStepUnit);// 2 <=> Cell Length
  streamTracer->SetInitialIntegrationStep(this->InitialIntegrationStep);//initial step length
  streamTracer->SetMinimumIntegrationStep(this->MinimumIntegrationStep);//Minimum Step Length
  streamTracer->SetMaximumIntegrationStep(this->MaximumIntegrationStep);//Maximum Step Length
  streamTracer->SetMaximumNumberOfSteps(this->MaximumNumberOfSteps);
  streamTracer->SetMaximumError(this->MaximumError);
  streamTracer->SetTerminalSpeed(this->TerminalSpeed);
  streamTracer->SetMaximumPropagation(dftMaxProp);
  streamTracer->SetSourceConnection(this->GetInputConnection(1,0));
  streamTracer->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, ArrayForGlyph); // idx==0  -> Vector selected
  streamTracer->Update();
  output->ShallowCopy(streamTracer->GetOutput());
  //
  vtkDataArray *arrToBeUsedToColor(output->GetPointData()->GetArray(ArrayForGlyph));
  vtkSmartPointer<vtkDataArray> arrColor(arrToBeUsedToColor->NewInstance());
  arrColor->ShallowCopy(arrToBeUsedToColor);
  arrColor->SetName(GetColorArrayName());
  int idx(output->GetPointData()->AddArray(arrColor));
  output->GetPointData()->SetActiveAttribute(idx,vtkDataSetAttributes::SCALARS);
  return 1;
}

const char vtkElectromagnetismStreamTraceur::NAME_COLOR_ARRAY[] = "Quantity To Display";

const char *vtkElectromagnetismStreamTraceur::GetColorArrayName()
{
  return NAME_COLOR_ARRAY;
}

//------------------------------------------------------------------------------
int vtkElectromagnetismStreamTraceur::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
  {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  }
  else if (port == 1)
  {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  }
  return 1;
}

//------------------------------------------------------------------------------
void vtkElectromagnetismStreamTraceur::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//------------------------------------------------------------------------------
vtkExecutive* vtkElectromagnetismStreamTraceur::CreateDefaultExecutive()
{
  return vtkCompositeDataPipeline::New();
}
