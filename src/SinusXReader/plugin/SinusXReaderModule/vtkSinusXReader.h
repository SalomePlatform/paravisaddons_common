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

#ifndef __vtkSinusXReader_h__
#define __vtkSinusXReader_h__

#include <vtkPolyDataAlgorithm.h>

class VTK_EXPORT vtkSinusXReader : public vtkPolyDataAlgorithm
{
public:
  static vtkSinusXReader *New();
  vtkTypeMacro(vtkSinusXReader, vtkPolyDataAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  vtkSetMacro(NumberOfSubdiv, int);
  vtkGetMacro(NumberOfSubdiv, int);

protected:
  vtkSinusXReader();
  ~vtkSinusXReader() override = default;

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  char *FileName = nullptr;
  int NumberOfSubdiv = 1;

private:
  vtkSinusXReader(const vtkSinusXReader &) = delete;
  void operator=(const vtkSinusXReader &) = delete;
};

#endif
