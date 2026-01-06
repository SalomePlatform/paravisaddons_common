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

#include "vtkElectromagnetismVecteur.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkMergeBlocks.h>
#include <vtkMultiBlockDataGroupFilter.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkPVGlyphFilter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCellCenters.h>
#include <vtkGlyphSource2D.h>

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkElectromagnetismVecteur);

const char vtkElectromagnetismVecteur::NAME_COLOR_ARRAY[] = "Quantity To Display";

const char *vtkElectromagnetismVecteur::GetColorArrayName()
{
  return NAME_COLOR_ARRAY;
}

 vtkElectromagnetismVecteur::vtkElectromagnetismVecteur():GlyphMode(ALL_POINTS)
  , MaximumNumberOfSamplePoints(5000)
  , Seed(1)
  , Stride(1)
  {
  }

// general_filters.xml : 1670
//-----------------------------------------------------------------------------
int vtkElectromagnetismVecteur::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  
  vtkDataArray *arr(this->GetInputArrayToProcess(1,inputVector));
  if(!arr)
  {
    vtkErrorMacro("No vector field selected in input !");
    return 0;
  }
  vtkInformation *inInfo( inputVector[0]->GetInformationObject(0) );
  vtkDataSet *input( vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT())) );
  const char *ArrayForGlyph(arr->GetName());
  vtkInformation *outInfo(outputVector->GetInformationObject(0));
  vtkPolyData *output(vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
  //
  vtkNew<vtkCellCenters> cc;
  cc->SetInputData(input);
  cc->Update();
  //
  vtkNew<vtkPVGlyphFilter> glyph;
  glyph->SetInputConnection(cc->GetOutputPort());
  glyph->SetGlyphMode(this->GlyphMode);       // vtkPVGlyphFilter::ALL_POINTS
  glyph->SetMaximumNumberOfSamplePoints(this->MaximumNumberOfSamplePoints);
  glyph->SetSeed(this->Seed);
  glyph->SetStride(this->Stride);
  glyph->SetVectorScaleMode(0); // vtkPVGlyphFilter::SCALE_BY_MAGNITUDE
  //
  vtkNew<vtkGlyphSource2D> arrow;
  arrow->SetGlyphTypeToArrow();
  arrow->SetFilled(false);
  glyph->SetSourceConnection(arrow->GetOutputPort());
  // idx,port,connection,fieldAssociation,name
  glyph->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, ""); // idx==0 -> scaleArray. "" means no scale array
  glyph->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, ArrayForGlyph); // idx==1 -> orientationArray
  glyph->SetScaleFactor(this->ScaleFactor);
  //
  glyph->Update();
  //
  output->ShallowCopy(glyph->GetOutput());
  //
  vtkDataArray *arrToBeUsedToColor(output->GetPointData()->GetArray(ArrayForGlyph));
  vtkSmartPointer<vtkDataArray> arrColor(arrToBeUsedToColor->NewInstance());
  arrColor->ShallowCopy(arrToBeUsedToColor);
  arrColor->SetName(GetColorArrayName());
  int idx(output->GetPointData()->AddArray(arrColor));
  output->GetPointData()->SetActiveAttribute(idx,vtkDataSetAttributes::SCALARS);
  return 1;
}

//-----------------------------------------------------------------------------
void vtkElectromagnetismVecteur::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
