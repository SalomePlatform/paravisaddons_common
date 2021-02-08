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
// Author : Anthony Geay (EDF R&D)

#ifndef vtkDepthVsTime_h__
#define vtkDepthVsTime_h__

#include <vtkDataObjectAlgorithm.h>

class vtkMutableDirectedGraph;

class VTK_EXPORT vtkDepthVsTime : public vtkDataObjectAlgorithm
{
public:
  static vtkDepthVsTime *New();
  vtkTypeMacro(vtkDepthVsTime, vtkDataObjectAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  void SetSourceData(vtkDataObject *input);

  void SetSourceConnection(vtkAlgorithmOutput *algOutput);

  vtkSetMacro(NbDiscrPtsAlongZ, int);
  vtkGetMacro(NbDiscrPtsAlongZ, int);

  int FillOutputPortInformation(int vtkNotUsed(port), vtkInformation *info) override;

protected:
  vtkDepthVsTime();
  ~vtkDepthVsTime();

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **inputVector, vtkInformationVector *vtkNotUsed(outputVector)) override;
  int RequestData(vtkInformation *request, vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int NbDiscrPtsAlongZ;
  int NumberOfTimeSteps;
  bool IsExecuting;
  int CurrentTimeIndex;
  class vtkInternal;
  vtkInternal *Internal;

private:
  vtkDepthVsTime(const vtkDepthVsTime &) = delete;
  void operator=(const vtkDepthVsTime &) = delete;
};

#endif
