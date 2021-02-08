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

#ifndef __vtkElectromagnetismVecteur_h__
#define __vtkElectromagnetismVecteur_h__

#include <vtkPolyDataAlgorithm.h>

#include <vtkSmartPointer.h>

#include <string>
#include <vector>

class vtkDoubleArray;
class vtkUnstructuredGrid;

class VTK_EXPORT vtkElectromagnetismVecteur : public vtkPolyDataAlgorithm
{
public:
  enum GlyphModeType
  {
    ALL_POINTS,
    EVERY_NTH_POINT,
    SPATIALLY_UNIFORM_DISTRIBUTION,
    SPATIALLY_UNIFORM_INVERSE_TRANSFORM_SAMPLING_SURFACE,
    SPATIALLY_UNIFORM_INVERSE_TRANSFORM_SAMPLING_VOLUME
  };
public:
  static vtkElectromagnetismVecteur* New();
  vtkTypeMacro(vtkElectromagnetismVecteur, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkGetMacro(ScaleFactor, double);
  vtkSetMacro(ScaleFactor, double);

  vtkSetClampMacro(GlyphMode, int, ALL_POINTS, SPATIALLY_UNIFORM_INVERSE_TRANSFORM_SAMPLING_VOLUME);
  vtkGetMacro(GlyphMode, int);

  vtkSetClampMacro(Stride, int, 1, VTK_INT_MAX);
  vtkGetMacro(Stride, int);

  vtkSetMacro(Seed, int);
  vtkGetMacro(Seed, int);

  vtkSetClampMacro(MaximumNumberOfSamplePoints, int, 1, VTK_INT_MAX);
  vtkGetMacro(MaximumNumberOfSamplePoints, int);

protected:
  vtkElectromagnetismVecteur();
  ~vtkElectromagnetismVecteur() override = default;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  static const char *GetColorArrayName();

  double ScaleFactor;

  int GlyphMode;
  int MaximumNumberOfSamplePoints;
  int Seed;
  int Stride;

private:
  vtkElectromagnetismVecteur(const vtkElectromagnetismVecteur&) = delete;
  void operator=(const vtkElectromagnetismVecteur&) = delete;

  static const char NAME_COLOR_ARRAY[];
};

#endif
