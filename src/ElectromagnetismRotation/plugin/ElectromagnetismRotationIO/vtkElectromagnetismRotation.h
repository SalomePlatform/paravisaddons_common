// Copyright (C) 2021-2024  CEA, EDF
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
// Author : Anthony Geay

#pragma once

#include "vtkMultiBlockDataSetAlgorithm.h"

#include <string>

class vtkMutableDirectedGraph;

class VTK_EXPORT vtkElectromagnetismRotation : public vtkMultiBlockDataSetAlgorithm
{
public:
  static vtkElectromagnetismRotation* New();
  vtkTypeMacro(vtkElectromagnetismRotation, vtkMultiBlockDataSetAlgorithm)
  void PrintSelf(ostream& os, vtkIndent indent);
  virtual int GetNumberOfGroupsFlagsArrays();
  const char *GetGroupsFlagsArrayName(int index);
  int GetGroupsFlagsArrayStatus(const char *name);
  virtual void SetGroupsFlagsStatus(const char *name, int status);
  void SetInsideOut(int val);
  void SetAxis(int axis);
  void SetAngularStep(char *angStep);

  // Description:
  // Every time the SIL is updated a this will return a different value.
  virtual int GetSILUpdateStamp();
  const char *GetMeshName();
  static const char* GetGrpStart();
  static const char* GetFamStart();
protected:
  vtkElectromagnetismRotation();
  ~vtkElectromagnetismRotation();

  int RequestInformation(vtkInformation *request,
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);

  int RequestData(vtkInformation *request, vtkInformationVector **inputVector,
      vtkInformationVector *outputVector);

  // Description:
  // This SIL stores the structure of the mesh/groups/cell types
  // that can be selected.
  virtual void SetSIL(vtkMutableDirectedGraph*);
  vtkGetObjectMacro(SIL, vtkMutableDirectedGraph);
protected:
  vtkMutableDirectedGraph *SIL;
  vtkTimeStamp SILTime;
private:
  vtkElectromagnetismRotation(const vtkElectromagnetismRotation&);
  void operator=(const vtkElectromagnetismRotation&); // Not implemented.
 private:
  //BTX
  //ETX
  class vtkElectromagnetismRotationInternal;
  vtkElectromagnetismRotationInternal *Internal;
  int InsideOut;
  int Axis;
  std::string AngularStep;
  mutable double RotationRotor;
};
