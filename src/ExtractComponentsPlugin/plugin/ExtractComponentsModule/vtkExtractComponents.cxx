// Copyright (C) 2021-2025  CEA, EDF
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
  Module:    vtkExtractComponents.cxx

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkExtractComponents.h"

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkFieldData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#ifdef WIN32
 #define NOMINMAX
 #include <algorithm>
#endif

vtkStandardNewMacro(vtkExtractComponents);

//----------------------------------------------------------------------------
vtkExtractComponents::vtkExtractComponents()
{
}

//----------------------------------------------------------------------------
vtkExtractComponents::~vtkExtractComponents()
{
  this->SetOutputArrayName(nullptr);
  this->ClearComponents();
}

//----------------------------------------------------------------------------
int vtkExtractComponents::FillInputPortInformation(int vtkNotUsed(port), vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

void vtkExtractComponents::SetGenerateVector(bool gvs)
{
  if (_gvs != gvs)
  {
    _gvs = gvs;
    this->Modified();
  }
}

//----------------------------------------------------------------------------
int vtkExtractComponents::RequestData(vtkInformation *vtkNotUsed(request),
                                      vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  // get the output info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->ShallowCopy(input);

  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);
  vtkInformation *info = this->GetInputArrayInformation(0);
  int fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_POINTS;
  if (info && info->Has(vtkDataObject::FIELD_ASSOCIATION()))
  {
    fieldAssociation = info->Get(vtkDataObject::FIELD_ASSOCIATION());
  }

  if (!inputArray)
  {
    vtkErrorMacro(<< "No data to extract");
    return 0;
  }

  if (this->InputArrayComponents.size() == 0)
  {
    return 1;
  }

  std::set<int>::const_iterator it = this->InputArrayComponents.begin();
  for (; it != this->InputArrayComponents.end(); ++it)
  {
    if (*it >= inputArray->GetNumberOfComponents() || *it < 0)
    {
      vtkErrorMacro(<< "Invalid component");
      return 0;
    }
  }

  if (!this->OutputArrayName)
  {
    vtkErrorMacro(<< "No output array name");
    return 0;
  }

  vtkSmartPointer<vtkDataArray> outputArray = inputArray->NewInstance();
  outputArray->SetName(this->OutputArrayName);
  int nbCompo(_gvs ? 3 : this->InputArrayComponents.size());
  outputArray->SetNumberOfComponents(nbCompo);
  outputArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());
  outputArray->CopyInformation(input->GetInformation());

  it = this->InputArrayComponents.begin();
  int imax(-1);
  for (int i = 0; it != this->InputArrayComponents.end(); ++it, i++)
  {
    if (!_gvs || i < 3)
    {
      imax = i;
      outputArray->CopyComponent(i, inputArray, *it);
      outputArray->SetComponentName(i, inputArray->GetComponentName(*it));
    }
  }
  if (_gvs)
  {
    for (int i = imax + 1; i < 3; i++)
    {
      vtkSmartPointer<vtkDataArray> tmpArray = inputArray->NewInstance();
      tmpArray->SetNumberOfComponents(1);
      tmpArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());
      tmpArray->Fill(0.);
      outputArray->CopyComponent(i, tmpArray, 0);
    }
    imax = 3;
  }
  imax = std::max(0, imax);
  if (fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_POINTS)
  {
    std::cerr << output->GetPointData()->GetNumberOfArrays() << std::endl;
    int idx(output->GetPointData()->AddArray(outputArray));
    output->GetPointData()->SetActiveAttribute(idx, imax > 1 ? vtkDataSetAttributes::VECTORS : vtkDataSetAttributes::SCALARS);
  }
  else if (fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_CELLS)
  {
    int idx(output->GetCellData()->AddArray(outputArray));
    output->GetCellData()->SetActiveAttribute(idx, imax > 1 ? vtkDataSetAttributes::VECTORS : vtkDataSetAttributes::SCALARS);
  }
  else if (fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_NONE)
  {
    output->GetFieldData()->AddArray(outputArray);
  }

  return 1;
}

//----------------------------------------------------------------------------
void vtkExtractComponents::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "InputArrayComponents: ";
  std::set<int>::iterator it;
  for (it = this->InputArrayComponents.begin(); it != this->InputArrayComponents.end(); ++it)
  {
    os << *it << " ";
  }
  os << std::endl;

  os << indent << "OutputArrayName: " << this->OutputArrayName << endl;
}

//----------------------------------------------------------------------------
void vtkExtractComponents::AddComponent(int component)
{
  this->InputArrayComponents.insert(component);
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkExtractComponents::ClearComponents()
{
  this->InputArrayComponents.clear();
  this->Modified();
}
