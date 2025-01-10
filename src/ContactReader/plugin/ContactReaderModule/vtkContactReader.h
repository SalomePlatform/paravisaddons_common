// Copyright (C) 2021-2025  CEA, EDF
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

#ifndef __vtkContactReader_h__
#define __vtkContactReader_h__

#include <vtkPolyDataAlgorithm.h>

class VTK_EXPORT vtkContactReader : public vtkPolyDataAlgorithm
{
public:
  static vtkContactReader *New();
  vtkTypeMacro(vtkContactReader, vtkPolyDataAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  void SetScaleFactor(double newScaleFactor);

protected:
  vtkContactReader();
  ~vtkContactReader() override;

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  char *FileName;
  double ScaleFactor;
  bool IsTimed;

private:
  vtkContactReader(const vtkContactReader &) = delete;
  void operator=(const vtkContactReader &) = delete;
};

#endif
