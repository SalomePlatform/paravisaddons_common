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

#ifndef vtkTemporalOnPoint_h__
#define vtkTemporalOnPoint_h__

#include <vtkDataObjectAlgorithm.h>

class vtkMutableDirectedGraph;

class VTK_EXPORT vtkTemporalOnPoint : public vtkDataObjectAlgorithm
{
public:
  static vtkTemporalOnPoint* New();
  vtkTypeMacro(vtkTemporalOnPoint, vtkDataObjectAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void SetSourceData(vtkDataObject* input);

  void SetSourceConnection(vtkAlgorithmOutput* algOutput);

  int FillOutputPortInformation(int, vtkInformation*) override;

protected:
  vtkTemporalOnPoint();
  ~vtkTemporalOnPoint() override;

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  int NumberOfTimeSteps;
  bool IsExecuting;
  int CurrentTimeIndex;
  class vtkInternal;
  vtkInternal* Internal;

private:
  vtkTemporalOnPoint(const vtkTemporalOnPoint&) = delete;
  void operator=(const vtkTemporalOnPoint&) = delete;
};

#endif
