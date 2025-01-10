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

#include "vtkGlyphCIH.h"

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

#include <cmath>
#include <sstream>

class MyGlyphException : public std::exception
{
public:
  MyGlyphException(const std::string& s):_reason(s) { }
  virtual const char *what() const throw() { return _reason.c_str(); }
  virtual ~MyGlyphException() throw() { }
private:
  std::string _reason;
};

//-----------------------------------------------------------------------------
void vtkGlyphCIH::ExtractInfo(
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
      vtkErrorMacro("Input dataSet must be a DataSet or single elt multi block dataset expected !");
      return;
    }
    if (input1->GetNumberOfBlocks() != 1)
    {
      vtkErrorMacro("Input dataSet is a multiblock dataset with not exactly one block ! Use "
                    "MergeBlocks or ExtractBlocks filter before calling this filter !");
      return;
    }
    vtkDataObject* input2(input1->GetBlock(0));
    if (!input2)
    {
      vtkErrorMacro("Input dataSet is a multiblock dataset with exactly one block but this single "
                    "element is NULL !");
      return;
    }
    vtkDataSet* input2c(vtkDataSet::SafeDownCast(input2));
    if (!input2c)
    {
      vtkErrorMacro(
        "Input dataSet is a multiblock dataset with exactly one block but this single element is "
        "not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
      return;
    }
    input = input2c;
  }

  if (!input)
  {
    vtkErrorMacro("Input data set is NULL !");
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
      usgIn = static_cast<vtkUnstructuredGrid *>(cd->GetOutput());
    }
    else
    {
      vtkNew<vtkMergeBlocks> filter;
      filter->SetMergePoints(0);
      filter->SetInputData(input1);
      filter->Update();
      usgIn = static_cast<vtkUnstructuredGrid *>(filter->GetOutput());
    }
  }
}

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkGlyphCIH);

int vtkGlyphCIH::FillInputPortInformation(int vtkNotUsed(port), vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

int vtkGlyphCIH::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

constexpr double radius_base = 0.03;
constexpr double radius_pointe = 0.1;
constexpr double z_base_pointe = 0.65;
constexpr double hauteur_pointe = 1.0 - z_base_pointe;
const double sqrt3Over2 = sqrt(3.) / 2;
const std::vector<double> fleche_base= {
  0.0,1.0,
  -sqrt3Over2,0.5,
  -sqrt3Over2,-0.5,
  0.0,-1.0,
  sqrt3Over2,-0.5,
  sqrt3Over2,0.5
  };
constexpr std::size_t N_CONN = 44;
constexpr std::size_t N_COORDS = 12;
constexpr std::size_t NB_CELL = 8;
constexpr std::size_t N_CONN_EFF = N_CONN-NB_CELL;
constexpr vtkIdType N_ACTIVE_CONN[N_CONN_EFF]={
  1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13,// 2 hexa
  15, 16, 17, 18, 20, 21, 22, 23, 25, 26, 27, 28, 30, 31, 32, 33, 35, 36, 37, 38, 40, 41, 42, 43// 6 quad4
};
constexpr vtkIdType N_CONN_STATIC[NB_CELL]={0, 7, 14, 19, 24, 29, 34, 39};
constexpr vtkIdType ref_conn[N_CONN]={6,0,1,2,3,4,5, 6,6,7,8,9,10,11, 4,0,1,7,6, 4,1,2,8,7, 4,2,3,9,8, 4,3,4,10,9, 4,4,5,11,10, 4,5,0,6,11};

double signOf(double val) { return val>=0?1.:-1.; }

template<class T>
void Algo(double widthFactor, vtkIdType inNbOfPts, double inNormMax, double scaleFactor, const double *inPts, double *coordsPtr, const T *inFieldDataPtr, vtkIdType *connPtr, T *fieldDataPtr)
{
  double sign(signOf(scaleFactor));
  scaleFactor = std::abs(scaleFactor);
  for(vtkIdType iPt = 0 ; iPt < inNbOfPts ; ++iPt, coordsPtr+=3*N_COORDS, connPtr+=N_CONN, fieldDataPtr+=3*NB_CELL)
    {
      for(auto i = 0 ; i < NB_CELL ; ++i)
        connPtr[N_CONN_STATIC[i]] = ref_conn[N_CONN_STATIC[i]];
      for(auto i = 0 ; i < N_CONN_EFF ; ++i)
        connPtr[N_ACTIVE_CONN[i]] = iPt*N_COORDS + ref_conn[N_ACTIVE_CONN[i]];
      T norm(vtkMath::Norm(inFieldDataPtr+3*iPt,3));
      //
      double mainZ[3] = {inFieldDataPtr[3*iPt+0]/norm,inFieldDataPtr[3*iPt+1]/norm,inFieldDataPtr[3*iPt+2]/norm};
      double mainX[3],mainY[3];
      vtkMath::Perpendiculars(mainZ,mainY,mainX,0.);
      double radiusBase(widthFactor*radius_base);
      double normRel(norm/inNormMax);
      // normRel entre 0. et 1.. 0 -> Pas de base de fleche. 1 fleche au taquet.
      /*
      double zPosBasePointe(normRel*scaleFactor*inNormMax-widthFactor*hauteur_pointe);
      double radiusPointe(widthFactor*radius_pointe);
      if(zPosBasePointe < 0.)
      {
        radiusPointe = (normRel*scaleFactor*inNormMax)/hauteur_pointe*radius_pointe;
        radiusBase = std::min(radiusBase,radiusPointe);
        zPosBasePointe = 0.0;
      }*/
      for(int i = 0 ; i < 6 ; ++i)
      {
        coordsPtr[3*i+0] = radiusBase*fleche_base[2*i+0]; coordsPtr[3*i+1] = radiusBase*fleche_base[2*i+1]; coordsPtr[3*i+2] = 0.0;
      }
      for(int i = 0 ; i < 6 ; ++i)
      {
        coordsPtr[6*3+3*i+0] = radiusBase*fleche_base[2*i+0]; coordsPtr[6*3+3*i+1] = radiusBase*fleche_base[2*i+1]; coordsPtr[6*3+3*i+2] = normRel*scaleFactor*inNormMax;//zPosBasePointe;
      }
      /*for(int i = 0 ; i < 6 ; ++i)
      {
        coordsPtr[12*3+3*i+0] = radiusPointe*fleche_base[2*i+0]; coordsPtr[12*3+3*i+1] = radiusPointe*fleche_base[2*i+1]; coordsPtr[12*3+3*i+2] = zPosBasePointe;
      }
      coordsPtr[18*3+0] = 0.; coordsPtr[18*3+1] = 0.; coordsPtr[18*3+2] = normRel*scaleFactor*inNormMax;*/
      // rotation
      for(auto i = 0 ; i < N_COORDS ; ++i )
        {
          double tmp[3] = {
            sign*(coordsPtr[i*3]*mainX[0]+coordsPtr[i*3+1]*mainY[0]+coordsPtr[i*3+2]*mainZ[0]),
            sign*(coordsPtr[i*3]*mainX[1]+coordsPtr[i*3+1]*mainY[1]+coordsPtr[i*3+2]*mainZ[1]),
            sign*(coordsPtr[i*3]*mainX[2]+coordsPtr[i*3+1]*mainY[2]+coordsPtr[i*3+2]*mainZ[2])
          };
          std::copy(tmp,tmp+3,coordsPtr+i*3);
        }
      // translation
      for(auto i = 0 ; i < N_COORDS ; ++i )
        { coordsPtr[i*3+0] += inPts[3*iPt+0]; coordsPtr[i*3+1] += inPts[3*iPt+1]; coordsPtr[i*3+2] += inPts[3*iPt+2]; }
      // field
      for(auto i = 0 ; i < NB_CELL ; ++i)
        { std::copy(inFieldDataPtr+3*iPt,inFieldDataPtr+3*(iPt+1),fieldDataPtr+i*3); }
    }
}

template<class T>
struct Traits
{
  using EltType = T;
};

template<>
struct Traits<double>
{
  using ArrayType = vtkDoubleArray;
};

template<>
struct Traits<float>
{
  using ArrayType = vtkFloatArray;
};


template<class T>
void Algo2(double widthFactor, double scaleFactor, vtkIdType inNbOfPts, const double *inPts, double *coordsPtr, vtkIdType *connPtr, typename Traits<T>::ArrayType *inFieldData, vtkDataArray *fieldData)
{
  
  const T *inFieldDataPtr(inFieldData->GetPointer(0));
  T *fieldDataPtr(reinterpret_cast<T *>(fieldData->GetVoidPointer(0)));
  //
  double inNormMin(std::numeric_limits<double>::max()),inNormMax(-std::numeric_limits<double>::max());
  for(vtkIdType i = 0 ; i < inNbOfPts ; ++i)
  {
    double norm(vtkMath::Norm(inFieldDataPtr+3*i,3));
    inNormMin = std::min(inNormMin,norm);
    inNormMax = std::max(inNormMax,norm);
  }
  //
  if(inNormMin == 0.)
    throw MyGlyphException("Norm min is null !");
  //
  Algo<T>(widthFactor,inNbOfPts,inNormMax,scaleFactor,inPts,coordsPtr,inFieldDataPtr,connPtr,fieldDataPtr);
}

void vtkGlyphCIH::SetInputArrayToProcess(int idx, int port, int connection, int ff, const char* name)
{
  if (idx == 0)
    this->FieldName = name;
  vtkPointSetAlgorithm::SetInputArrayToProcess(idx, port, connection, ff, name);
}

//-----------------------------------------------------------------------------
int vtkGlyphCIH::RequestData(vtkInformation* vtkNotUsed(request),vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  try
  {
    vtkInformation* outInfo(outputVector->GetInformationObject(0));
    vtkPointSet* output(vtkPointSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
    vtkSmartPointer<vtkUnstructuredGrid> usgIn;
    this->ExtractInfo(inputVector[0], usgIn);
    //
    vtkNew<vtkCellCenters> cc;
    cc->SetInputData(usgIn);
    cc->Update();
    vtkDataArray *inFieldDataGen(cc->GetOutput()->GetPointData()->GetArray(this->FieldName.c_str()));
    vtkIdType inNbOfPts(cc->GetOutput()->GetNumberOfPoints());
    const double *inPts(static_cast<double *>(cc->GetOutput()->GetPoints()->GetData()->GetVoidPointer(0)));
    //
    if(!inFieldDataGen)
    {
      std::ostringstream oss; oss << "No such point field with name \"" << this->FieldName << "\" !";
      throw MyGlyphException(oss.str());
    }
    vtkSmartPointer<vtkDataArray> fieldData;
    fieldData.TakeReference(inFieldDataGen->NewInstance());
    fieldData->SetName(this->FieldName.c_str());
    for(int i=0;i<3;i++)
      fieldData->SetComponentName(i,inFieldDataGen->GetComponentName(i));
    fieldData->SetNumberOfComponents(3);
    fieldData->SetNumberOfTuples(inNbOfPts*NB_CELL);
    //
    vtkNew<vtkDoubleArray> coords;
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(inNbOfPts*N_COORDS);
    vtkNew<vtkIdTypeArray> conn;
    conn->SetNumberOfComponents(1);
    conn->SetNumberOfTuples(inNbOfPts*N_CONN);
    double *coordsPtr(coords->GetPointer(0));
    vtkIdType *connPtr(conn->GetPointer(0));
    //
    /*double inMaxCellSize(0);
    {
      vtkNew<vtkMeshQuality> mq;
      mq->SetTriangleQualityMeasureToArea();
      mq->SetQuadQualityMeasureToArea();
      mq->SetInputData(usgIn);
      mq->Update();
      vtkDoubleArray *area(vtkDoubleArray::SafeDownCast(mq->GetOutput()->GetCellData()->GetArray("Quality")));
      inMaxCellSize = area->GetMaxNorm();
    }*/
    //
    if(vtkDoubleArray::SafeDownCast(inFieldDataGen))
    {
      vtkDoubleArray *inFieldData = vtkDoubleArray::SafeDownCast(inFieldDataGen);
      Algo2<double>(this->WidthFactor,this->ScaleFactor,inNbOfPts,inPts,coordsPtr,connPtr,inFieldData,fieldData);
    }
    else if(vtkFloatArray::SafeDownCast(inFieldDataGen))
    {
      vtkFloatArray *inFieldData = vtkFloatArray::SafeDownCast(inFieldDataGen);
      Algo2<float>(this->WidthFactor,this->ScaleFactor,inNbOfPts,inPts,coordsPtr,connPtr,inFieldData,fieldData);
    }
    else
    {
      throw MyGlyphException("Only float64 and float32 managed for input point field !");
    }
    //
    vtkNew<vtkPolyData> pd;
    vtkNew<vtkCellArray> cb;
    cb->SetCells(inNbOfPts*NB_CELL,conn);
    pd->SetPolys(cb);
    //
    vtkNew<vtkPoints> pts;
    pts->SetData(coords);
    pd->SetPoints(pts);
    // on ajoute tous les pointsarrays
    vtkCellData *pdc(pd->GetCellData());
    pdc->AddArray(fieldData);
    //
    output->ShallowCopy(pd);
  }
  catch(MyGlyphException& e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkComplexMode::RequestInformation : " << e.what() << std::endl;
    if(this->HasObserver("ErrorEvent") )
      this->InvokeEvent("ErrorEvent",const_cast<char *>(oss.str().c_str()));
    else
      vtkOutputWindowDisplayErrorText(const_cast<char *>(oss.str().c_str()));
    vtkObject::BreakOnError();
  }
  return 1;
}

//-----------------------------------------------------------------------------
void vtkGlyphCIH::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
