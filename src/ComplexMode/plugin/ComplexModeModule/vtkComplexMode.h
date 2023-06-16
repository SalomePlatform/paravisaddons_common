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
// Author : Anthony Geay (EDF R&D)

#ifndef vtkComplexMode_h__
#define vtkComplexMode_h__

#include <vtkUnstructuredGridAlgorithm.h>

class vtkMutableDirectedGraph;

class VTK_EXPORT vtkComplexMode : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkComplexMode* New();
  vtkTypeMacro(vtkComplexMode, vtkUnstructuredGridAlgorithm)
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkGetMacro(Factor,double);
  vtkSetClampMacro(Factor,double,0.,VTK_DOUBLE_MAX);

  vtkGetMacro(Phase,double);
  vtkSetClampMacro(Phase,double,-180.,180.);

  vtkGetMacro(AnimationTime,double);
  vtkSetClampMacro(AnimationTime,double,0.,1.);

  void SetInputArrayToProcess(int idx, int port, int connection, int fieldAssociation, const char *name);

protected:
  vtkComplexMode();
  ~vtkComplexMode();

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector, vtkInformationVector *outputVector);

  int RequestData(vtkInformation *request, vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector);

private:
  vtkComplexMode(const vtkComplexMode&);
  void operator=(const vtkComplexMode&); // Not implemented.

protected:
  double Factor;
  double Phase;
  double AnimationTime;
  class vtkComplexModeInternal;
  vtkComplexModeInternal* Internal;
};

#endif
