// Copyright (C) 2021-2026  CEA, EDF
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

#ifndef vtkSphereAlongLines_h__
#define vtkSphereAlongLines_h__

#include <vtkPolyDataAlgorithm.h>

class VTK_EXPORT vtkSphereAlongLines : public vtkPolyDataAlgorithm
{
public:
  static vtkSphereAlongLines* New();
  vtkTypeMacro(vtkSphereAlongLines, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkGetMacro(AnimationTime, double);
  vtkSetClampMacro(AnimationTime, double, 0., 1.);

  vtkGetMacro(WalkType, int);
  vtkSetMacro(WalkType, int);

protected:
  vtkSphereAlongLines();
  ~vtkSphereAlongLines() override;

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  class vtkInternals;
  vtkInternals* Internal;
  double AnimationTime;
  int WalkType;

private:
  vtkSphereAlongLines(const vtkSphereAlongLines&) = delete;
  void operator=(const vtkSphereAlongLines&) = delete;
};

#endif
