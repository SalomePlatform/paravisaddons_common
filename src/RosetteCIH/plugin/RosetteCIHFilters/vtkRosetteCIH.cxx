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

#include "vtkRosetteCIH.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkMergeBlocks.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDoubleArray.h>
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
#include <vtkWarpVector.h>

#include <stdexcept>
#include <sstream>

void RosetteCIHEmitThrowInternal(const std::string& msg)
{
  //vtkErrorMacro(msg);
  throw std::runtime_error( msg );
}

#define RosetteCIHEmitThrow(text)                  \
{                                                  \
    vtkErrorMacro(text);                           \
    std::ostringstream oss; oss << text;           \
    RosetteCIHEmitThrowInternal(oss.str());        \
}

//-----------------------------------------------------------------------------
void vtkRosetteCIH::ExtractInfo(
  vtkInformationVector* inputVector, vtkSmartPointer<vtkUnstructuredGrid>& usgIn)
{
  vtkInformation* inputInfo(inputVector->GetInformationObject(0));
  vtkDataSet* input(0);
  vtkDataSet* input0(vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  vtkMultiBlockDataSet* input1(
    vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  if (input0)
  {
    input = input0;
  }
  else
  {
    if (!input1)
    {
      RosetteCIHEmitThrow("Input dataSet must be a DataSet or single elt multi block dataset expected !");
      return;
    }
    if (input1->GetNumberOfBlocks() != 1)
    {
      RosetteCIHEmitThrow("Input dataSet is a multiblock dataset with not exactly one block ! Use "
                    "MergeBlocks or ExtractBlocks filter before calling this filter !");
      return;
    }
    vtkDataObject* input2(input1->GetBlock(0));
    if (!input2)
    {
      RosetteCIHEmitThrow("Input dataSet is a multiblock dataset with exactly one block but this single "
                    "element is NULL !");
      return;
    }
    vtkDataSet* input2c(vtkDataSet::SafeDownCast(input2));
    if (!input2c)
    {
      RosetteCIHEmitThrow(
        "Input dataSet is a multiblock dataset with exactly one block but this single element is "
        "not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
      return;
    }
    input = input2c;
  }

  if (!input)
  {
    RosetteCIHEmitThrow("Input data set is NULL !");
    return;
  }

  usgIn = vtkUnstructuredGrid::SafeDownCast(input);
  if (!usgIn)
  {
    if (!input1)
    {
      vtkNew<vtkMultiBlockDataGroupFilter> mb;
      vtkNew<vtkMergeBlocks> cd;
      mb->AddInputData(input);
      cd->SetInputConnection(mb->GetOutputPort());
      cd->SetMergePoints(0);
      cd->Update();
      usgIn = static_cast<vtkUnstructuredGrid*>(cd->GetOutput());
    }
    else
    {
      vtkNew<vtkMergeBlocks> filter;
      filter->SetMergePoints(0);
      filter->SetInputData(input1);
      filter->Update();
      usgIn = static_cast<vtkUnstructuredGrid*>(filter->GetOutput());
    }
  }
}

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkRosetteCIH);

//-----------------------------------------------------------------------------
vtkSmartPointer<vtkDataSet> vtkRosetteCIH::GenerateGlyphLinesFor(vtkUnstructuredGrid* usgIn,
  const char* keyPoint, const char COMPRESS_TRACTION[])
{
  vtkFieldData* dsa(usgIn->GetCellData());
  std::string arrayForGlyph(this->RetrieveFieldForGlyph(usgIn, keyPoint));
  int compoId(-1);
  vtkDoubleArray* arrayForPosNeg2(RetrieveFieldForPost(usgIn, keyPoint, compoId));
  // vtkAbstractArray *arrayForPosNeg(dsa->GetAbstractArray(FIELD_NAME_2));
  // vtkDoubleArray *arrayForPosNeg2(vtkDoubleArray::SafeDownCast(arrayForPosNeg));
  int nbCompo(arrayForPosNeg2->GetNumberOfComponents());
  vtkIdType nbTuples(arrayForPosNeg2->GetNumberOfTuples());
  vtkNew<vtkDoubleArray> compressionOrTraction;
  compressionOrTraction->SetNumberOfComponents(1);
  compressionOrTraction->SetNumberOfTuples(nbTuples);
  compressionOrTraction->SetName(COMPRESS_TRACTION);
  const double* pt(arrayForPosNeg2->GetPointer(0));
  double* ptOut(compressionOrTraction->GetPointer(0));
  for (vtkIdType i = 0; i < nbTuples; i++)
  {
    if (pt[i * nbCompo + compoId] > 0.)
      ptOut[i] = 1.;
    else
      ptOut[i] = -1.;
  }
  int arrId(dsa->AddArray(compressionOrTraction));
  //
  vtkNew<vtkPVGlyphFilter> glyph;
  glyph->SetInputData(usgIn);
  glyph->SetGlyphMode(0);       // vtkPVGlyphFilter::ALL_POINTS
  glyph->SetVectorScaleMode(0); // vtkPVGlyphFilter::SCALE_BY_MAGNITUDE

  //
  vtkNew<vtkLineSource> arrow;
  glyph->SetSourceConnection(arrow->GetOutputPort());
  // idx,port,connection,fieldAssociation,name
  glyph->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arrayForGlyph.c_str()); // idx==0 -> scaleArray
  glyph->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS,
    arrayForGlyph.c_str()); // idx==1 -> orientationArray
  glyph->SetScaleFactor(this->ScaleFactor);
  glyph->Update();

  return vtkSmartPointer<vtkDataSet>(glyph->GetOutput());
}

//-----------------------------------------------------------------------------
void vtkRosetteCIH::PostTraitementT1etT2(
  vtkUnstructuredGrid* usgIn, vtkUnstructuredGrid* output)
{
  constexpr char COMPRESS_TRACTION[] = "CompressionOrTraction";
  // "RESUNL__SIRO_ELEM_T1_Vector" , "RESUNL__SIRO_ELEM_T1"
  vtkSmartPointer<vtkDataSet> gl1 =
    this->GenerateGlyphLinesFor(usgIn, "T1", COMPRESS_TRACTION);
  vtkSmartPointer<vtkDataSet> gl2 =
    this->GenerateGlyphLinesFor(usgIn, "T2", COMPRESS_TRACTION);
  //
  vtkNew<vtkDataSetSurfaceFilter> surface;
  surface->SetNonlinearSubdivisionLevel(0);
  surface->SetInputData(usgIn);
  surface->Update();
  vtkNew<vtkPolyData> surfaceCpy;
  surfaceCpy->ShallowCopy(surface->GetOutput());
  vtkNew<vtkDoubleArray> compressionOrTraction;
  auto nbOfTuples(surface->GetOutput()->GetNumberOfPoints());
  compressionOrTraction->SetNumberOfComponents(1);
  compressionOrTraction->SetNumberOfTuples(nbOfTuples);
  compressionOrTraction->SetName(COMPRESS_TRACTION);
  compressionOrTraction->Fill(NAN);
  surfaceCpy->GetPointData()->AddArray(compressionOrTraction);
  //
  vtkNew<vtkMultiBlockDataGroupFilter> mb;
  vtkNew<vtkMergeBlocks> cd;
  mb->AddInputData(surfaceCpy);
  mb->AddInputData(gl1);
  mb->AddInputData(gl2);
  cd->SetInputConnection(mb->GetOutputPort());
  cd->SetMergePoints(0);
  cd->Update();
  //
  output->ShallowCopy(cd->GetOutput());
  //
  vtkFieldData* dsa(output->GetPointData());
  int nbOfArrays(dsa->GetNumberOfArrays());
  for (int i = nbOfArrays - 1; i >= 0; i--)
  {
    const char* arrName(dsa->GetArrayName(i));
    if (std::string(arrName) != COMPRESS_TRACTION)
    {
      dsa->RemoveArray(i);
    }
  }
  output->GetPointData()->SetActiveAttribute(0, vtkDataSetAttributes::SCALARS);
}

//-----------------------------------------------------------------------------
int vtkRosetteCIH::ComponentIdOfArray(vtkAbstractArray* array, const std::string& compoName)
{
  int nbCompo(array->GetNumberOfComponents());
  int ret(-1);
  for (int i = 0; i < nbCompo; i++)
  {
    if (compoName == array->GetComponentName(i))
    {
      if (ret != -1)
      {
        RosetteCIHEmitThrow("ComponentIdOfArray : already found !");
        return ret;
      }
      ret = i;
    }
  }
  if (ret == -1)
  {
    RosetteCIHEmitThrow(
      "ComponentIdOfArray : component " << compoName << " in array " << array->GetName() << " !");
  }
  return ret;
}

//-----------------------------------------------------------------------------
std::string vtkRosetteCIH::GenerateAValidFieldForGlyph(
  vtkUnstructuredGrid* usgInCpy, const std::string& arrayName, const char* keyPoint)
{
  vtkFieldData* dsa(usgInCpy->GetCellData());
  vtkAbstractArray* array(dsa->GetAbstractArray(arrayName.c_str()));
  vtkDoubleArray* array2(vtkDoubleArray::SafeDownCast(array));
  //
  std::string ret(arrayName);
  ret += std::string("_vveeccttoorr");
  int compoIds[3] = { -1, -1, -1 };
  for (int i = 0; i < 3; i++)
  {
    std::string compoName("SIG_");
    compoName += keyPoint;
    compoName += 'X' + i;
    compoIds[i] = ComponentIdOfArray(array, compoName);
  }
  //
  vtkIdType nbTuples(array2->GetNumberOfTuples());
  int nbCompo(array2->GetNumberOfComponents());
  vtkNew<vtkDoubleArray> vect;
  vect->SetNumberOfComponents(3);
  vect->SetNumberOfTuples(nbTuples);
  vect->SetName(ret.c_str());

  const double* pt(array2->GetPointer(0));
  double* ptOut(vect->GetPointer(0));
  for (vtkIdType i = 0; i < nbTuples; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      ptOut[3 * i + j] = pt[nbCompo * i + compoIds[j]];
    }
  }
  //
  dsa->AddArray(vect);
  return ret;
}

//-----------------------------------------------------------------------------
bool vtkRosetteCIH::EndWith(const std::string& arrayName, const std::string& end)
{
  std::size_t lenOfLastChance(end.length());
  if (arrayName.length() < lenOfLastChance)
  {
    return false;
  }

  std::string endOfArrayName(arrayName.substr(arrayName.length() - lenOfLastChance));
  return endOfArrayName == end;
}

//-----------------------------------------------------------------------------
bool vtkRosetteCIH::IsFirstChance(const std::string& arrayName, const char* keyPoint)
{
  std::string PATTERN("SIRO_ELEM");
  PATTERN += std::string("_") + keyPoint;
  return this->EndWith(arrayName, PATTERN);
}

//-----------------------------------------------------------------------------
bool vtkRosetteCIH::IsLastChanceArray(const std::string& arrayName)
{
  return this->EndWith(arrayName, "SIRO_ELEM");
}

//-----------------------------------------------------------------------------
std::string vtkRosetteCIH::GetFieldName(vtkUnstructuredGrid* usgInCpy, const char* keyPoint)
{
  vtkFieldData* dsa(usgInCpy->GetCellData());
  std::string arrayNameOK;
  int nbOfArrays(dsa->GetNumberOfArrays());
  bool found(false);
  for (int i = 0; i < nbOfArrays; ++i)
  {
    vtkAbstractArray* arrayAbstract(dsa->GetAbstractArray(i));
    std::string arrayName(arrayAbstract->GetName());
    if (this->IsFirstChance(arrayName, keyPoint) || this->IsLastChanceArray(arrayName))
    {
      if (found)
      {
        RosetteCIHEmitThrow("GetFieldName : already found !");
      }
      arrayNameOK = arrayName;
      found = true;
    }
  }
  if (!found)
  {
    RosetteCIHEmitThrow("GetFieldName : Impossible to find a valid array !");
  }
  return arrayNameOK;
}

//-----------------------------------------------------------------------------
std::string vtkRosetteCIH::RetrieveFieldForGlyph(
  vtkUnstructuredGrid* usgInCpy, const char* keyPoint)
{
  std::string arrayNameOK(this->GetFieldName(usgInCpy, keyPoint));
  return this->GenerateAValidFieldForGlyph(usgInCpy, arrayNameOK, keyPoint);
}

//-----------------------------------------------------------------------------
vtkDoubleArray* vtkRosetteCIH::RetrieveFieldForPost(
  vtkUnstructuredGrid* usgInCpy, const char* keyPoint, int& compId)
{
  std::string FIELD_NAME_2(this->GetFieldName(usgInCpy, keyPoint));
  vtkFieldData* dsa(usgInCpy->GetCellData());
  vtkAbstractArray* arrayForPosNeg(dsa->GetAbstractArray(FIELD_NAME_2.c_str()));
  vtkDoubleArray* arrayForPosNeg2(vtkDoubleArray::SafeDownCast(arrayForPosNeg));
  std::string compoToFind("SIG_");
  compoToFind += keyPoint;
  compId = this->ComponentIdOfArray(arrayForPosNeg, compoToFind);
  return arrayForPosNeg2;
}

//-----------------------------------------------------------------------------
void vtkRosetteCIH::PostTraitementOnlyOneCompo(vtkUnstructuredGrid* usgIn,
  vtkUnstructuredGrid* output, const char* keyPoint,
  const char* COMPRESS_TRACTION)
{
  vtkNew<vtkUnstructuredGrid> usgInCpy;
  usgInCpy->DeepCopy(usgIn);
  //
  vtkFieldData* dsa(usgInCpy->GetCellData());
  int compId(-1);
  // vtkAbstractArray *arrayForPosNeg(dsa->GetAbstractArray(FIELD_NAME_2));
  std::string arrayForGlyph(this->RetrieveFieldForGlyph(usgInCpy, keyPoint));
  vtkDoubleArray* arrayForPosNeg2(this->RetrieveFieldForPost(usgInCpy, keyPoint, compId));
  // vtkDoubleArray::SafeDownCast(arrayForPosNeg);
  int nbCompo(arrayForPosNeg2->GetNumberOfComponents());
  vtkIdType nbTuples(arrayForPosNeg2->GetNumberOfTuples());

  vtkNew<vtkDoubleArray> compressionOrTraction;
  compressionOrTraction->SetNumberOfComponents(1);
  compressionOrTraction->SetNumberOfTuples(nbTuples);
  compressionOrTraction->SetName(COMPRESS_TRACTION);

  const double* pt(arrayForPosNeg2->GetPointer(0));
  double* ptOut(compressionOrTraction->GetPointer(0));
  double valMin(std::numeric_limits<double>::max());
  double valMax(-std::numeric_limits<double>::max());

  for (vtkIdType i = 0; i < nbTuples; i++)
  {
    double val(pt[i * nbCompo + compId]);
    valMin = std::min(valMin, val);
    valMax = std::max(valMax, val);
    ptOut[i] = val;
  }
  //
  for (int i = dsa->GetNumberOfArrays() - 1; i >= 0; i--)
  {
    if (arrayForGlyph != dsa->GetAbstractArray(i)->GetName())
    {
      dsa->RemoveArray(i);
    }
  }
  int arrId(dsa->AddArray(compressionOrTraction));

  vtkNew<vtkLineSource> arrow;

  vtkNew<vtkDataSetSurfaceFilter> surface;
  surface->SetNonlinearSubdivisionLevel(0);
  surface->SetInputData(usgInCpy);

  vtkNew<vtkPolyDataNormals> normals;
  normals->ComputeCellNormalsOn();
  normals->ComputePointNormalsOn();
  normals->SplittingOff();
  normals->SetInputConnection(surface->GetOutputPort());
  normals->Update();

  // for some reasons, the glyph filter removes scalars and normals, we have to duplicate them
  vtkDataArray* normalsArray = normals->GetOutput()->GetCellData()->GetNormals();

  vtkSmartPointer<vtkDataArray> savedNormalsArray;
  savedNormalsArray.TakeReference(normalsArray->NewInstance());
  savedNormalsArray->DeepCopy(normalsArray);
  savedNormalsArray->SetName("CellNormals");

  normals->GetOutput()->GetCellData()->AddArray(savedNormalsArray);

  vtkNew<vtkWarpVector> warp;
  warp->SetInputConnection(normals->GetOutputPort());
  warp->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Normals");
  warp->SetScaleFactor(normals->GetOutput()->GetLength()/1000);
  warp->Update();

  vtkNew<vtkPVGlyphFilter> glyph;
  glyph->SetInputConnection(warp->GetOutputPort());
  glyph->SetGlyphMode(0);       // vtkPVGlyphFilter::ALL_POINTS
  glyph->SetVectorScaleMode(0); // vtkPVGlyphFilter::SCALE_BY_MAGNITUDE
  glyph->SetSourceConnection(arrow->GetOutputPort());
  // idx,port,connection,fieldAssociation,name
  glyph->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arrayForGlyph.c_str()); // idx==0 -> scaleArray
  glyph->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS,
    arrayForGlyph.c_str()); // idx==1 -> orientationArray
  glyph->SetScaleFactor(this->ScaleFactor);

  vtkNew<vtkRibbonFilter> ribbon;
  ribbon->SetWidth(this->WidthFactor);
  ribbon->VaryWidthOff();
  ribbon->UseDefaultNormalOff();
  ribbon->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "CellNormals");
  ribbon->SetInputConnection(glyph->GetOutputPort());
  ribbon->Update();

  vtkDataSet* ret = ribbon->GetOutput();

  vtkFieldData* fieldData = ret->GetPointData();
  for (int i = fieldData->GetNumberOfArrays() - 1; i >= 0; i--)
  {
    fieldData->RemoveArray(i);
  }

  fieldData = ret->GetCellData();
  for (int i = fieldData->GetNumberOfArrays() - 1; i >= 0; i--)
  {
    fieldData->RemoveArray(i);
  }

  vtkNew<vtkDoubleArray> compressionOrTractionNaN;
  compressionOrTractionNaN->SetNumberOfComponents(1);
  compressionOrTractionNaN->SetNumberOfTuples(ret->GetNumberOfCells());
  compressionOrTractionNaN->SetName(COMPRESS_TRACTION);
  compressionOrTractionNaN->Fill(NAN);
  fieldData->AddArray(compressionOrTractionNaN);

  vtkNew<vtkMultiBlockDataGroupFilter> mb;
  mb->AddInputData(usgInCpy);
  mb->AddInputData(ret);

  vtkNew<vtkMergeBlocks> cd;
  cd->SetInputConnection(mb->GetOutputPort());
  cd->SetMergePoints(0);
  cd->Update();

  output->ShallowCopy(cd->GetOutput());

  int arrayId;
  output->GetCellData()->GetAbstractArray(COMPRESS_TRACTION, arrayId);
  output->GetCellData()->SetActiveAttribute(arrayId, vtkDataSetAttributes::SCALARS);
}

//-----------------------------------------------------------------------------
void vtkRosetteCIH::PostTraitementT1(vtkUnstructuredGrid* usgIn, vtkUnstructuredGrid* output)
{
  // constexpr char FIELD_NAME[]="RESUNL__SIRO_ELEM_T1_Vector";
  // constexpr char FIELD_NAME_2[]="RESUNL__SIRO_ELEM_T1";
  constexpr char COMPRESS_TRACTION[] = "Contrainte specifique 1";
  this->PostTraitementOnlyOneCompo(usgIn, output, "T1", COMPRESS_TRACTION);
}

//-----------------------------------------------------------------------------
void vtkRosetteCIH::PostTraitementT2(vtkUnstructuredGrid* usgIn, vtkUnstructuredGrid* output)
{
  // constexpr char FIELD_NAME[]="RESUNL__SIRO_ELEM_T2_Vector";
  // constexpr char FIELD_NAME_2[]="RESUNL__SIRO_ELEM_T2";
  constexpr char COMPRESS_TRACTION[] = "Contrainte specifique 3";
  this->PostTraitementOnlyOneCompo(usgIn, output, "T2", COMPRESS_TRACTION);
}

//-----------------------------------------------------------------------------
int vtkRosetteCIH::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkInformation* outInfo(outputVector->GetInformationObject(0));
  vtkUnstructuredGrid* output(
    vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
  //
  vtkSmartPointer<vtkUnstructuredGrid> usgIn;
  try
  {
    this->ExtractInfo(inputVector[0], usgIn);
    switch (this->TypeOfDisplay)
    {
      case 0:
        this->PostTraitementT1etT2(usgIn, output);
        break;
      case 1:
        this->PostTraitementT1(usgIn, output);
        break;
      case 2:
        this->PostTraitementT2(usgIn, output);
        break;
      default:
        RosetteCIHEmitThrow("GetFieldName : Impossible to find a valid array !");
    }
  }
  catch(const std::exception& e)
  {
    std::cerr << "vtkRosetteCIH::RequestData : " << e.what() << std::endl;
    return 0;
  }

  return 1;
}

//-----------------------------------------------------------------------------
void vtkRosetteCIH::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
