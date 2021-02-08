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
  Module:    vtkSMMyNumberOfComponentsDomain.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkSMMyNumberOfComponentsDomain
 * @brief   int range domain based on the number of
 * components available in a particular data array.
 *
 * vtkSMMyNumberOfComponentsDomain is used for properties that allow the user to
 * choose the component number (or associated name) to process for the choosen array.
 * It needs two required properties with following functions:
 * * Input -- input property for the filter.
 * * ArraySelection -- string vector property used to select the array.
 * This domain will not work if either of the required properties is missing.
*/

#ifndef vtkSMMyNumberOfComponentsDomain_h
#define vtkSMMyNumberOfComponentsDomain_h

#include <vtkSMEnumerationDomain.h>

class vtkSMSourceProxy;
class vtkSMInputArrayDomain;

class VTK_EXPORT vtkSMMyNumberOfComponentsDomain : public vtkSMEnumerationDomain
{
public:
  static vtkSMMyNumberOfComponentsDomain* New();
  vtkTypeMacro(vtkSMMyNumberOfComponentsDomain, vtkSMEnumerationDomain);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Updates the range based on the scalar range of the currently selected
   * array. This requires Input (vtkSMProxyProperty) and ArraySelection
   * (vtkSMStringVectorProperty) properties. Currently, this uses
   * only the first component of the array.
   */
  virtual void Update(vtkSMProperty* prop);

protected:
  vtkSMMyNumberOfComponentsDomain();
  ~vtkSMMyNumberOfComponentsDomain() override;

  /**
   * Internal update method doing the actual work.
   */
  void Update(
    const char* arrayname, int arrayType, vtkSMSourceProxy* sp, vtkSMInputArrayDomain* iad, int outputport);

private:
  vtkSMMyNumberOfComponentsDomain(const vtkSMMyNumberOfComponentsDomain&) = delete;
  void operator=(const vtkSMMyNumberOfComponentsDomain&) = delete;
};

#endif
