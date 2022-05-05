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
// Author : Anthony Geay (EDF R&D)

#pragma once

#include "vtkDataObjectAlgorithm.h"

class vtkPolyData;

class VTK_EXPORT vtkSpatialPfl : public vtkDataObjectAlgorithm
{
public:
  static vtkSpatialPfl* New();
  vtkTypeMacro(vtkSpatialPfl, vtkDataObjectAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void SetSourceData(vtkDataObject* input);

  void SetSourceConnection(vtkAlgorithmOutput* algOutput);

  //@{
  /**
   * Set/Get the number of points that will be considered along the polyline source
   * before the probing.
   */
  vtkGetMacro(ResampleInput, bool);
  vtkSetMacro(ResampleInput, bool);
  vtkBooleanMacro(ResampleInput, bool);
  //@}

  //@{
  /**
   * Set/Get the number of points that will be considered along the polyline source
   * before the probing.
   */
  vtkGetMacro(NumberOfSamples, int);
  vtkSetMacro(NumberOfSamples, int);
  //@}

  int FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info) override;

protected:
  vtkSpatialPfl();
  ~vtkSpatialPfl() override = default;

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  vtkPolyData* ResamplePolyLine(vtkPolyData* pl);

  bool ResampleInput = false;
  int NumberOfSamples = 0;

private:
  vtkSpatialPfl(const vtkSpatialPfl&) = delete;
  void operator=(const vtkSpatialPfl&) = delete;
};
