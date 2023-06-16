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

/*=========================================================================

  Program:   ParaView
  Module:    vtkExtractComponents.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkExtractComponents
 * @brief   Extract a component of an attribute.
 *
 * vtkExtractComponents Extract a component of an attribute.
*/

#ifndef vtkExtractComponents_h
#define vtkExtractComponents_h

#include <vtkDataSetAlgorithm.h>

#include <set>

class VTK_EXPORT vtkExtractComponents : public vtkDataSetAlgorithm
{
public:
  static vtkExtractComponents *New();
  vtkTypeMacro(vtkExtractComponents, vtkDataSetAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  void AddComponent(int);
  void ClearComponents();

  vtkSetStringMacro(OutputArrayName);
  vtkGetStringMacro(OutputArrayName);

  void SetGenerateVector(bool gvs);

protected:
  vtkExtractComponents();
  ~vtkExtractComponents() override;

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  virtual int FillInputPortInformation(int port, vtkInformation *info) override;

  std::set<int> InputArrayComponents;
  char *OutputArrayName = nullptr;
  bool _gvs = false;

private:
  vtkExtractComponents(const vtkExtractComponents &) = delete;
  void operator=(const vtkExtractComponents &) = delete;
};

#endif
