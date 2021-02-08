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

  Program:   ParaView
  Module:    vtkSMMyNumberOfComponentsDomain.cxx

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSMMyNumberOfComponentsDomain.h"

#include <vtkObjectFactory.h>
#include <vtkPVArrayInformation.h>
#include <vtkPVDataInformation.h>
#include <vtkPVDataSetAttributesInformation.h>
#include <vtkSMDomainIterator.h>
#include <vtkSMInputArrayDomain.h>
#include <vtkSMInputProperty.h>
#include <vtkSMSourceProxy.h>
#include <vtkSMStringVectorProperty.h>

#include <sstream>

vtkStandardNewMacro(vtkSMMyNumberOfComponentsDomain);
//----------------------------------------------------------------------------
vtkSMMyNumberOfComponentsDomain::vtkSMMyNumberOfComponentsDomain()
{
}

//----------------------------------------------------------------------------
vtkSMMyNumberOfComponentsDomain::~vtkSMMyNumberOfComponentsDomain()
{
}

//----------------------------------------------------------------------------
void vtkSMMyNumberOfComponentsDomain::Update(vtkSMProperty*)
{
  vtkSMProxyProperty* ip = vtkSMProxyProperty::SafeDownCast(this->GetRequiredProperty("Input"));
  vtkSMStringVectorProperty* svp =
    vtkSMStringVectorProperty::SafeDownCast(this->GetRequiredProperty("ArraySelection"));
  if (!ip || !svp)
  {
    // Missing required properties.
    this->RemoveAllEntries();
    this->DomainModified();
    return;
  }

 /* if (svp->GetNumberOfUncheckedElements() != 5 && svp->GetNumberOfUncheckedElements() != 2 &&
    svp->GetNumberOfUncheckedElements() != 1)
  {
    // We can only handle array selection properties with 5, 2 or 1 elements.
    // For 5 elements the array name is at indices [4]; for 2
    // elements it's at [1], while for 1 elements, it's at [0].
    this->RemoveAllEntries();
    return;
  }*/

  int index = svp->GetNumberOfUncheckedElements() - 1;
  const char* arrayName = svp->GetUncheckedElement(index);
  int arrayType = atoi(svp->GetUncheckedElement(index - 1));
  if (!arrayName || arrayName[0] == 0)
  {
    // No array choosen.
    this->RemoveAllEntries();
    this->DomainModified();
    return;
  }

  vtkSMInputArrayDomain* iad = nullptr;
  vtkSMDomainIterator* di = ip->NewDomainIterator();
  di->Begin();
  while (!di->IsAtEnd())
  {
    // We have to figure out whether we are working with cell data,
    // point data or both.
    iad = vtkSMInputArrayDomain::SafeDownCast(di->GetDomain());
    if (iad)
    {
      break;
    }
    di->Next();
  }
  di->Delete();
  if (!iad)
  {
    // Failed to locate a vtkSMInputArrayDomain on the input property, which is
    // required.
    this->RemoveAllEntries();
    this->DomainModified();
    return;
  }

  vtkSMInputProperty* inputProp = vtkSMInputProperty::SafeDownCast(ip);
  unsigned int numProxs = ip->GetNumberOfUncheckedProxies();
  for (unsigned int i = 0; i < numProxs; i++)
  {
    // Use the first input
    vtkSMSourceProxy* source = vtkSMSourceProxy::SafeDownCast(ip->GetUncheckedProxy(i));
    if (source)
    {
      this->Update(arrayName, arrayType, source, iad,
        (inputProp ? inputProp->GetUncheckedOutputPortForConnection(i) : 0));
      return;
    }
  }
}

//---------------------------------------------------------------------------
void vtkSMMyNumberOfComponentsDomain::Update(
  const char* arrayName, int arrayType, vtkSMSourceProxy* sp,
  vtkSMInputArrayDomain* iad, int outputport)
{
  // Make sure the outputs are created.
  sp->CreateOutputPorts();
  this->RemoveAllEntries();
  vtkPVDataInformation* info = sp->GetDataInformation(outputport);
  if (!info)
  {
    this->DomainModified();
    return;
  }

  vtkPVArrayInformation* ai = 0;
  if (arrayType == vtkDataObject::FIELD_ASSOCIATION_POINTS)
  {
    ai = info->GetPointDataInformation()->GetArrayInformation(arrayName);
  }
  else if (arrayType == vtkDataObject::FIELD_ASSOCIATION_CELLS)
  {
    ai = info->GetCellDataInformation()->GetArrayInformation(arrayName);
  }
  else if (arrayType == vtkDataObject::FIELD_ASSOCIATION_VERTICES)
  {
    ai = info->GetVertexDataInformation()->GetArrayInformation(arrayName);
  }
  else if (arrayType == vtkDataObject::FIELD_ASSOCIATION_EDGES)
  {
    ai = info->GetEdgeDataInformation()->GetArrayInformation(arrayName);
  }
  else if (arrayType == vtkDataObject::FIELD_ASSOCIATION_ROWS)
  {
    ai = info->GetRowDataInformation()->GetArrayInformation(arrayName);
  }
  else if (arrayType == vtkDataObject::FIELD_ASSOCIATION_NONE)
  {
    ai = info->GetFieldDataInformation()->GetArrayInformation(arrayName);
  }

  if (ai)
  {
    for (int i = 0; i < ai->GetNumberOfComponents(); i++)
    {
      const char* name = ai->GetComponentName(i);
      std::ostringstream cname;
      if (!name || name[0] == 0)
      {
        cname << i;
      }
      else
      {
        cname << name;
      }
      this->AddEntry(cname.str().c_str(), i);
    }
  }
  this->DomainModified();
}

//----------------------------------------------------------------------------
void vtkSMMyNumberOfComponentsDomain::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
