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
// Author : Anthony Geay (EDF R&D)

#pragma once

#include <vtkUnstructuredGridAlgorithm.h>
//#include <vtkPolyDataAlgorithm.h>

class vtkMutableDirectedGraph;
class vtkImplicitFunction;
class vtkCylinder;

class VTK_EXPORT vtkElectromagnetismFluxDisc : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkElectromagnetismFluxDisc* New();
  vtkTypeMacro(vtkElectromagnetismFluxDisc, vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void SetCutFunction(vtkImplicitFunction* func);

  vtkMTimeType GetMTime();

  vtkSetClampMacro(RadialResolution, int, 1, VTK_INT_MAX);
  vtkGetMacro(RadialResolution, int);
  
  vtkSetClampMacro(CircumferentialResolution, int, 3, VTK_INT_MAX);
  vtkGetMacro(CircumferentialResolution, int);

protected:
  vtkElectromagnetismFluxDisc();
  ~vtkElectromagnetismFluxDisc();
  

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  vtkCylinder* _cyl;
  class vtkInternals;
  vtkInternals* Internal;

private:
  vtkElectromagnetismFluxDisc(const vtkElectromagnetismFluxDisc&) = delete;
  void operator=(const vtkElectromagnetismFluxDisc&) = delete;
  
  int RadialResolution;
  int CircumferentialResolution;
};
