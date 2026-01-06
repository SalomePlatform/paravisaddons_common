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

#ifndef __vtkExplodePolyLine_h__
#define __vtkExplodePolyLine_h__

#include <vtkPointSetAlgorithm.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkMergeBlocks.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLineSource.h>
#include <vtkMultiBlockDataGroupFilter.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkRibbonFilter.h>
#include <vtkPVGlyphFilter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkTable.h>
#include <vtkTessellatorFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVariant.h>
#include <vtkVariantArray.h>
#include <vtkMeshQuality.h>
#include <vtkCellCenters.h>

#include <string>
#include <vector>

class vtkDoubleArray;

class VTK_EXPORT vtkExplodePolyLine : public vtkPointSetAlgorithm
{
public:
  static vtkExplodePolyLine* New();
  vtkTypeMacro(vtkExplodePolyLine, vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  int FillInputPortInformation(int vtkNotUsed(port), vtkInformation *info) override;
  int FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info) override;

protected:
  vtkExplodePolyLine() = default;
  ~vtkExplodePolyLine() override = default;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  void ExtractInfo(vtkInformationVector* inputVector, vtkSmartPointer<vtkPolyData>& usgIn);

private:
  vtkExplodePolyLine(const vtkExplodePolyLine&) = delete;
  void operator=(const vtkExplodePolyLine&) = delete;
};

#endif
