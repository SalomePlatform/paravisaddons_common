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

#ifndef __vtkGlyphCIH_h__
#define __vtkGlyphCIH_h__

#include <vtkPointSetAlgorithm.h>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <string>
#include <vector>

class vtkDoubleArray;

class VTK_EXPORT vtkGlyphCIH : public vtkPointSetAlgorithm
{
public:
  static vtkGlyphCIH* New();
  vtkTypeMacro(vtkGlyphCIH, vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkGetMacro(ScaleFactor, double);
  vtkSetMacro(ScaleFactor, double);

  vtkGetMacro(WidthFactor, double);
  vtkSetMacro(WidthFactor, double);

  int FillInputPortInformation(int vtkNotUsed(port), vtkInformation *info) override;
  int FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info) override;
  
  void SetInputArrayToProcess(int idx, int port, int connection, int fieldAssociation, const char* name) override;

protected:
  vtkGlyphCIH() = default;
  ~vtkGlyphCIH() override = default;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  void ExtractInfo(vtkInformationVector* inputVector, vtkSmartPointer<vtkUnstructuredGrid>& usgIn);

  double ScaleFactor;
  double WidthFactor;
  int TypeOfDisplay;
  std::string FieldName;

private:
  vtkGlyphCIH(const vtkGlyphCIH&) = delete;
  void operator=(const vtkGlyphCIH&) = delete;
};

#endif
