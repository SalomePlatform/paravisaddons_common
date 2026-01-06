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
// Author : Anthony Geay

#ifndef vtkExtractThreeD_h__
#define vtkExtractThreeD_h__

#include <vtkDataSetAlgorithm.h>

class vtkMutableDirectedGraph;

class VTK_EXPORT vtkExtractThreeD : public vtkDataSetAlgorithm
{
public:
  static vtkExtractThreeD *New();
  vtkTypeMacro(vtkExtractThreeD, vtkDataSetAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent) override;

protected:
  vtkExtractThreeD();
  ~vtkExtractThreeD();

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

  int RequestData(vtkInformation *request, vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  class vtkExtractThreeDInternal;
  vtkExtractThreeDInternal *Internal;

private:
  vtkExtractThreeD(const vtkExtractThreeD &) = delete;
  void operator=(const vtkExtractThreeD &) = delete;
};

#endif
