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

#include "vtkTorseurCIH.h"

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAlgorithmOutput.h>
#include <vtkCompositeDataToUnstructuredGridFilter.h>
#include <vtkDataArraySelection.h>
#include <vtkDataObjectTreeIterator.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationDataObjectKey.h>
#include <vtkInformationDataObjectMetaDataKey.h>
#include <vtkInformationDoubleVectorKey.h>
#include <vtkInformationQuadratureSchemeDefinitionVectorKey.h>
#include <vtkInformationStringKey.h>
#include <vtkInformationVector.h>
#include <vtkLongArray.h>
#include <vtkMultiBlockDataGroupFilter.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkObjectFactory.h>
#include <vtkQuadratureSchemeDefinition.h>
#include <vtkStringArray.h>
#include <vtkTable.h>
#include <vtkUnsignedCharArray.h>

#include <deque>
#include <map>
#include <set>
#include <sstream>

vtkStandardNewMacro(vtkTorseurCIH);
///////////////////

void ExtractInfo(vtkInformationVector* inputVector, vtkSmartPointer<vtkUnstructuredGrid>& usgIn)
{
  vtkInformation* inputInfo(inputVector->GetInformationObject(0));
  vtkDataSet* input = nullptr;
  vtkDataSet* input0(vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  vtkMultiBlockDataSet* input1(
    vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  if (input0)
    input = input0;
  else
  {
    if (!input1)
      return ;
    if (input1->GetNumberOfBlocks() != 1)
      return ;
    vtkDataObject* input2(input1->GetBlock(0));
    if (!input2)
      return ;
    vtkDataSet* input2c(vtkDataSet::SafeDownCast(input2));
    if (!input2c)
      return ;
    input = input2c;
  }
  if (!input)
    return ;
  usgIn.TakeReference(vtkUnstructuredGrid::SafeDownCast(input));
  if (!usgIn.Get())
  {
    if (!input1)
    {
      vtkNew<vtkMultiBlockDataGroupFilter> mb;
      vtkNew<vtkCompositeDataToUnstructuredGridFilter> cd;
      mb->AddInputData(input);
      cd->SetInputConnection(mb->GetOutputPort());
      cd->SetMergePoints(0);
      cd->Update();
      usgIn = cd->GetOutput();
    }
    else
    {
      vtkNew<vtkCompositeDataToUnstructuredGridFilter> filter;
      filter->SetMergePoints(0);
      filter->SetInputData(input1);
      filter->Update();
      vtkUnstructuredGrid* res(filter->GetOutput());
      usgIn.TakeReference(res);
      if (res)
        res->Register(nullptr);
    }
  }
  else
    usgIn->Register(nullptr);
}


vtkSmartPointer<vtkTable> ComputeTorseurCIH(vtkUnstructuredGrid* usgIn)
{
  double area = (double) usgIn->GetNumberOfCells();
  double InertiaNormale = (double) usgIn->GetNumberOfPoints();
  double inertia = 2 * area;
  double inertiaOther = 3 * InertiaNormale;
  double centerOfMass[3] = {0.,0.,0.};
  double ForceNormale[3] = {1.,1.,1.};
  double TangentForce[3] = {2.,2.,2.};
  double normalFace[3] = {3.,3.,3.};
  double outputAxis[3] = {4.,4.,4.};
  double tangentOther[3] = {5.,5.,5.};
  double momentum[3] = {6.,6.,6.};
  //
  vtkSmartPointer<vtkTable> ret(vtkSmartPointer<vtkTable>::New());
  vtkSmartPointer<vtkStringArray> col0(vtkSmartPointer<vtkStringArray>::New());
  constexpr int NB_ROWS = 11;
  col0->SetNumberOfComponents(1);
  col0->SetNumberOfTuples(NB_ROWS);
  col0->SetName("Grandeur");
  // scalaire
  col0->SetValue(0, strdup("Aire"));
  col0->SetValue(1, strdup("Inertie Normal"));
  col0->SetValue(2, strdup("Inertie Tangentielle principale"));
  col0->SetValue(3, strdup("Inertie Tangentielle secondaire"));
  // vectoriel
  col0->SetValue(4, strdup("Position du centre de gravite"));
  col0->SetValue(5, strdup("Effort Normal"));
  col0->SetValue(6, strdup("Effort Tangentiel"));
  col0->SetValue(7, strdup("Axe Normal"));
  col0->SetValue(8, strdup("Axe Tangentiel principal"));
  col0->SetValue(9, strdup("Axe Tangentiel secondaire"));
  col0->SetValue(10, strdup("Moment au centre de gravite"));
  ret->AddColumn(col0);
  //
  vtkSmartPointer<vtkDoubleArray> col1(vtkSmartPointer<vtkDoubleArray>::New());
  col1->SetName("X");
  col1->SetNumberOfComponents(1);
  col1->SetNumberOfTuples(NB_ROWS);
  col1->SetValue(0, area);
  col1->SetValue(1, InertiaNormale);
  col1->SetValue(2, inertia);
  col1->SetValue(3, inertiaOther);
  col1->SetValue(4, centerOfMass[0]);
  col1->SetValue(5, ForceNormale[0]);
  col1->SetValue(6, TangentForce[0]);
  col1->SetValue(7, normalFace[0]);
  col1->SetValue(8, outputAxis[0]);
  col1->SetValue(9, tangentOther[0]);
  col1->SetValue(10, momentum[0]);
  ret->AddColumn(col1);
  //
  vtkSmartPointer<vtkDoubleArray> col2(vtkSmartPointer<vtkDoubleArray>::New());
  col2->SetName("Y");
  col2->SetNumberOfComponents(1);
  col2->SetNumberOfTuples(NB_ROWS);
  col2->SetValue(0, 0.);
  col2->SetValue(1, 0.);
  col2->SetValue(2, 0.);
  col2->SetValue(3, 0.);
  col2->SetValue(4, centerOfMass[1]);
  col2->SetValue(5, ForceNormale[1]);
  col2->SetValue(6, TangentForce[1]);
  col2->SetValue(7, normalFace[1]);
  col2->SetValue(8, outputAxis[1]);
  col2->SetValue(9, tangentOther[1]);
  col2->SetValue(10, momentum[1]);
  ret->AddColumn(col2);
  //
  vtkSmartPointer<vtkDoubleArray> col3(vtkSmartPointer<vtkDoubleArray>::New());
  col3->SetName("Z");
  col3->SetNumberOfComponents(1);
  col3->SetNumberOfTuples(NB_ROWS);
  col3->SetValue(0, 0.);
  col3->SetValue(1, 0.);
  col3->SetValue(2, 0.);
  col3->SetValue(3, 0.);
  col3->SetValue(4, centerOfMass[2]);
  col3->SetValue(5, ForceNormale[2]);
  col3->SetValue(6, TangentForce[2]);
  col3->SetValue(7, normalFace[2]);
  col3->SetValue(8, outputAxis[2]);
  col3->SetValue(9, tangentOther[2]);
  col3->SetValue(10, momentum[2]);
  ret->AddColumn(col3);
  return ret;
}

////////////////////

int vtkTorseurCIH::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
  return 1;
}

int vtkTorseurCIH::RequestInformation(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkSmartPointer<vtkUnstructuredGrid> usgIn;
  return 1;
}

int vtkTorseurCIH::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  {
    vtkSmartPointer<vtkUnstructuredGrid> usgIn;
    ExtractInfo(inputVector[0], usgIn);
    //
    vtkSmartPointer<vtkTable> ret(ComputeTorseurCIH(usgIn));
    vtkInformation* inInfo(inputVector[0]->GetInformationObject(0));
    vtkInformation* outInfo(outputVector->GetInformationObject(0));
    vtkTable* output(vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
    output->ShallowCopy(ret);
  }
  return 1;
}

void vtkTorseurCIH::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
