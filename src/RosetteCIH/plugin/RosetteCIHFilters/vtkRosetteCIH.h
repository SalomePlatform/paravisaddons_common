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

#ifndef __vtkRosetteCIH_h__
#define __vtkRosetteCIH_h__

#include <vtkUnstructuredGridAlgorithm.h>

#include <vtkSmartPointer.h>

#include <string>
#include <vector>

class vtkDoubleArray;

class VTK_EXPORT vtkRosetteCIH : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkRosetteCIH* New();
  vtkTypeMacro(vtkRosetteCIH, vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkGetMacro(ScaleFactor, double);
  vtkSetMacro(ScaleFactor, double);

  vtkGetMacro(WidthFactor, double);
  vtkSetMacro(WidthFactor, double);

  vtkGetMacro(TypeOfDisplay, int);
  vtkSetMacro(TypeOfDisplay, int);

protected:
  vtkRosetteCIH() = default;
  ~vtkRosetteCIH() override = default;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  void ExtractInfo(vtkInformationVector* inputVector, vtkSmartPointer<vtkUnstructuredGrid>& usgIn);

  vtkSmartPointer<vtkDataSet> GenerateGlyphLinesFor(
    vtkUnstructuredGrid* usgIn, const char* keyPoint, const char COMPRESS_TRACTION[]);

  void PostTraitementT1etT2(vtkUnstructuredGrid* usgIn, vtkUnstructuredGrid* output);

  int ComponentIdOfArray(vtkAbstractArray* array, const std::string& compoName);

  std::string GenerateAValidFieldForGlyph(
    vtkUnstructuredGrid* usgInCpy, const std::string& arrayName, const char* keyPoint);

  bool EndWith(const std::string& arrayName, const std::string& end);

  bool IsFirstChance(const std::string& arrayName, const char* keyPoint);
  bool IsLastChanceArray(const std::string& arrayName);

  std::string GetFieldName(vtkUnstructuredGrid* usgInCpy, const char* keyPoint);

  std::string RetrieveFieldForGlyph(vtkUnstructuredGrid* usgInCpy, const char* keyPoint);

  vtkDoubleArray* RetrieveFieldForPost(
    vtkUnstructuredGrid* usgInCpy, const char* keyPoint, int& compId);

  void PostTraitementOnlyOneCompo(vtkUnstructuredGrid* usgIn, vtkUnstructuredGrid* output,
    const char* keyPoint, const char* COMPRESS_TRACTION);

  void PostTraitementT1(vtkUnstructuredGrid* usgIn, vtkUnstructuredGrid* output);
  void PostTraitementT2(vtkUnstructuredGrid* usgIn, vtkUnstructuredGrid* output);

  double ScaleFactor;
  double WidthFactor;
  int TypeOfDisplay;

private:
  vtkRosetteCIH(const vtkRosetteCIH&) = delete;
  void operator=(const vtkRosetteCIH&) = delete;
};

#endif
