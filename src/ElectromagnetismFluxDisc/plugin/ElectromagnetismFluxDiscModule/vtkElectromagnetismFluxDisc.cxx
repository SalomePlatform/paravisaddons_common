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
// Author : Anthony Geay (EDF R&D)

#include "vtkElectromagnetismFluxDisc.h"

#include "vtkAdjacentVertexIterator.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkCylinder.h"
#include "vtkNew.h"
#include "vtkCutter.h"
#include "vtkTransform.h"

#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"
#include  "vtkMultiBlockDataSet.h"

#include "vtkInformationStringKey.h"
#include "vtkAlgorithmOutput.h"
#include "vtkObjectFactory.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkDataSet.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataArraySelection.h"
#include "vtkTimeStamp.h"
#include "vtkInEdgeIterator.h"
#include "vtkInformationDataObjectKey.h"
#include "vtkExecutive.h"
#include "vtkVariantArray.h"
#include "vtkStringArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkDemandDrivenPipeline.h"
#include "vtkDataObjectTreeIterator.h"
#include "vtkWarpScalar.h"
#include "vtkDiskSource.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkResampleWithDataSet.h"
#include "vtkPointDataToCellData.h"
#include "vtkMeshQuality.h"

#include <cmath>


#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include <map>
#include <deque>
#include <sstream>
#include <algorithm>

vtkStandardNewMacro(vtkElectromagnetismFluxDisc);

class VTK_EXPORT MZCException : public std::exception
{
public:
  MZCException(const std::string& s):_reason(s) { }
  virtual const char *what() const throw() { return _reason.c_str(); }
  virtual ~MZCException() throw() { }
private:
  std::string _reason;
};


void ExtractInfo(vtkInformationVector *inputVector, vtkDataSet *& usgIn)
{
  vtkInformation *inputInfo(inputVector->GetInformationObject(0));
  vtkDataSet *input(0);
  vtkDataSet *input0(vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  vtkMultiBlockDataSet *input1(vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  if(input0)
    input=input0;
  else
    {
      if(!input1)
        throw MZCException("Input dataSet must be a DataSet or single elt multi block dataset expected !");
      if(input1->GetNumberOfBlocks()!=1)
      {
        std::cerr << "**** " << input1->GetNumberOfBlocks() << std::endl;
        throw MZCException("Input dataSet is a multiblock dataset with not exactly one block ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
      }
      vtkDataObject *input2(input1->GetBlock(0));
      if(!input2)
        throw MZCException("Input dataSet is a multiblock dataset with exactly one block but this single element is NULL !");
      vtkDataSet *input2c(vtkDataSet::SafeDownCast(input2));
      if(!input2c)
        throw MZCException("Input dataSet is a multiblock dataset with exactly one block but this single element is not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
      input=input2c;
    }
  if(!input)
    throw MZCException("Input data set is NULL !");
  vtkPointData *att(input->GetPointData());
  if(!att)
    throw MZCException("Input dataset has no point data attribute ! Impossible to deduce a developed surface on it !");
  usgIn=input;
}

class vtkElectromagnetismFluxDisc::vtkInternals
{
public:
  vtkNew<vtkCutter> Cutter;
};

////////////////////

vtkElectromagnetismFluxDisc::vtkElectromagnetismFluxDisc():_cyl(nullptr),Internal(new vtkInternals),RadialResolution(80),CircumferentialResolution(80)
{
}

vtkElectromagnetismFluxDisc::~vtkElectromagnetismFluxDisc()
{
  delete this->Internal;
}

int vtkElectromagnetismFluxDisc::RequestInformation(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{ 
  //std::cerr << "########################################## vtkElectromagnetismFluxDisc::RequestInformation ##########################################" << std::endl;
  try
    {
      vtkDataSet *usgIn(0);
      ExtractInfo(inputVector[0],usgIn);
    }
  catch(MZCException& e)
    {
      vtkErrorMacro("Exception has been thrown in vtkElectromagnetismFluxDisc::RequestInformation : " << e.what());
      return 0;
    }
  return 1;
}

double ComputeFlux(vtkIdType nbOfTuples, const double *area, const double *vector3Field, const double axis[3])
{
  double ret(0.0);
  for( vtkIdType i = 0 ; i < nbOfTuples ; ++i )
    ret += area[i] * vtkMath::Dot(vector3Field+i*3,axis);
  return ret;
}

template<class T>
void Rotate3DAlg(const double *center, const double *vect, double angle, vtkIdType nbNodes, const T *coordsIn, T *coordsOut)
{
  double sina(sin(angle));
  double cosa(cos(angle));
  double vectorNorm[3];
  T matrix[9];
  T matrixTmp[9];
  double norm(sqrt(vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]));
  if(norm<std::numeric_limits<double>::min())
    throw MZCException("Rotate3DAlg : magnitude of input vector is too close of 0. !");
  std::transform(vect,vect+3,vectorNorm,std::bind2nd(std::multiplies<double>(),1/norm));
  //rotation matrix computation
  matrix[0]=cosa; matrix[1]=0.; matrix[2]=0.; matrix[3]=0.; matrix[4]=cosa; matrix[5]=0.; matrix[6]=0.; matrix[7]=0.; matrix[8]=cosa;
  matrixTmp[0]=vectorNorm[0]*vectorNorm[0]; matrixTmp[1]=vectorNorm[0]*vectorNorm[1]; matrixTmp[2]=vectorNorm[0]*vectorNorm[2];
  matrixTmp[3]=vectorNorm[1]*vectorNorm[0]; matrixTmp[4]=vectorNorm[1]*vectorNorm[1]; matrixTmp[5]=vectorNorm[1]*vectorNorm[2];
  matrixTmp[6]=vectorNorm[2]*vectorNorm[0]; matrixTmp[7]=vectorNorm[2]*vectorNorm[1]; matrixTmp[8]=vectorNorm[2]*vectorNorm[2];
  std::transform(matrixTmp,matrixTmp+9,matrixTmp,std::bind2nd(std::multiplies<T>(),1-cosa));
  std::transform(matrix,matrix+9,matrixTmp,matrix,std::plus<T>());
  matrixTmp[0]=0.; matrixTmp[1]=-vectorNorm[2]; matrixTmp[2]=vectorNorm[1];
  matrixTmp[3]=vectorNorm[2]; matrixTmp[4]=0.; matrixTmp[5]=-vectorNorm[0];
  matrixTmp[6]=-vectorNorm[1]; matrixTmp[7]=vectorNorm[0]; matrixTmp[8]=0.;
  std::transform(matrixTmp,matrixTmp+9,matrixTmp,std::bind2nd(std::multiplies<T>(),sina));
  std::transform(matrix,matrix+9,matrixTmp,matrix,std::plus<T>());
  //rotation matrix computed.
  T tmp[3];
  for(vtkIdType i=0; i<nbNodes; i++)
    {
      std::transform(coordsIn+i*3,coordsIn+(i+1)*3,center,tmp,std::minus<T>());
      coordsOut[i*3]=matrix[0]*tmp[0]+matrix[1]*tmp[1]+matrix[2]*tmp[2]+(T)center[0];
      coordsOut[i*3+1]=matrix[3]*tmp[0]+matrix[4]*tmp[1]+matrix[5]*tmp[2]+(T)center[1];
      coordsOut[i*3+2]=matrix[6]*tmp[0]+matrix[7]*tmp[1]+matrix[8]*tmp[2]+(T)center[2];
    }
}

int vtkElectromagnetismFluxDisc::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkElectromagnetismFluxDisc::RequestData        ##########################################" << std::endl;
  try
    {
      if(!_cyl)
        throw MZCException("No cylinder object as cut function !");
      double center[3],axis[3],radius,orthoAxis[3];
      const double ZVec[3] = {0.,0.,1.};
      vtkAbstractTransform* trf(_cyl->GetTransform());
      {
        _cyl->GetCenter(center);
        _cyl->GetAxis(axis[0],axis[1],axis[2]);
        radius=_cyl->GetRadius();
      }
      if(trf)
        {
          double axis3[3]={center[0]+0.,center[1]+1.,center[2]+0.},axis4[3];
          trf->TransformPoint(axis3,axis4);
          std::transform(axis4,axis4+3,center,axis,[](double a, double b) { return b-a; });
          axis[1]=-axis[1];
          if(std::isnan(axis[0]) && std::isnan(axis[1]) && std::isnan(axis[2]))
            { axis[0]=0.; axis[1]=-1.; axis[2]=0.; }
        }
      
      //std::cerr << trf << " jjj " << axis[0] << " " << axis[1] << " " << axis[2]  << " : " << center[0] <<  " " << center[1] << " " << center[2] <<  " " " " << " -> " << radius << std::endl;
      vtkDataSet *usgIn(0);
      ExtractInfo(inputVector[0],usgIn);
      //
      vtkNew<vtkUnstructuredGrid> outputMesh;
      {
        vtkIdType nbPoints(this->RadialResolution*this->CircumferentialResolution+1);
        vtkNew<vtkDoubleArray> coords;
        coords->SetNumberOfComponents(3); coords->SetNumberOfTuples(nbPoints);
        double *coordsPtr(coords->GetPointer(0));
        coordsPtr[0] = 0; coordsPtr[1] = 0; coordsPtr[2] = 0; coordsPtr+=3;
        for(int circI = 0 ; circI < this->CircumferentialResolution ; ++circI)
        {
          double a(2*M_PI*(double(circI)/double(this->CircumferentialResolution)));
          double c(cos(a)),s(sin(a));
          for(int radI = 0 ; radI < this->RadialResolution ; ++radI)
          {
            coordsPtr[0] = c*(double(radI+1)/double(this->RadialResolution))*radius;
            coordsPtr[1] = s*(double(radI+1)/double(this->RadialResolution))*radius;
            coordsPtr[2] = 0;
            coordsPtr+=3;
          }
        }
        vtkNew<vtkPoints> pts;
        pts->SetData(coords);
        outputMesh->SetPoints(pts);
        //
        vtkIdType nbOfCells(this->CircumferentialResolution*this->RadialResolution);
        vtkNew<vtkCellArray> cells;
        vtkNew<vtkIdTypeArray> cellsData;
        vtkNew<vtkIdTypeArray> cellLocations;
        vtkNew<vtkUnsignedCharArray> cellTypes;
        //
        cellTypes->SetNumberOfComponents(1); cellTypes->SetNumberOfTuples(nbOfCells);
        cellLocations->SetNumberOfComponents(1); cellLocations->SetNumberOfTuples(nbOfCells);
        cellsData->SetNumberOfComponents(1); cellsData->SetNumberOfTuples( ( 5*(this->RadialResolution)-1 ) * this->CircumferentialResolution );
        vtkIdType *clPtr(cellLocations->GetPointer(0)),*cdPtr(cellsData->GetPointer(0));
        vtkIdType offset(0),deltaPt(this->RadialResolution);
        unsigned char *ctPtr(cellTypes->GetPointer(0));
        for( int iCirc = 0 ; iCirc < this->CircumferentialResolution - 1; ++iCirc )
        {
          vtkIdType zeDelta(iCirc*deltaPt);
          for( int iRadial = 0 ; iRadial < this->RadialResolution ; ++iRadial )
          {
            *clPtr++ = offset;
            if(iRadial!=0)
            {
              *ctPtr++ = VTK_QUAD;
              cdPtr[0] = 4 ; cdPtr[1] = zeDelta + iRadial ; cdPtr[2] = zeDelta + iRadial+1;
              cdPtr[3] = (zeDelta + deltaPt + iRadial+1); cdPtr[4] = ( zeDelta + deltaPt + iRadial);
              cdPtr+=5;
              offset += 4;
            }
            else
            {
              *ctPtr++ = VTK_TRIANGLE;
              cdPtr[0] = 3 ; cdPtr[1] = 0 ; cdPtr[2] = zeDelta + 1; cdPtr[3] = (zeDelta + deltaPt +1);
              cdPtr += 4;
              offset += 3;
            }
          }
        }
        vtkIdType zeDelta((this->CircumferentialResolution - 1)*deltaPt);
        for( int iRadial = 0 ; iRadial < this->RadialResolution ; ++iRadial )
          {
            *clPtr++ = offset;
            if(iRadial!=0)
            {
              *ctPtr++ = VTK_QUAD;
              cdPtr[0] = 4 ; cdPtr[1] = zeDelta + iRadial ; cdPtr[2] = zeDelta + iRadial+1;
              cdPtr[3] = iRadial+1; cdPtr[4] = iRadial;
              cdPtr+=5;
              offset += 4;
            }
            else
            {
              *ctPtr++ = VTK_TRIANGLE;
              cdPtr[0] = 3 ; cdPtr[1] = 0 ; cdPtr[2] = zeDelta + 1; cdPtr[3] = 1;
              cdPtr += 4;
              offset += 3;
            }
          }
        //
        cells->SetCells(nbOfCells,cellsData);
        outputMesh->SetCells(cellTypes,cellLocations,cells);
      }
      // Rotation
      {
        vtkIdType nbPoints(outputMesh->GetNumberOfPoints());
        vtkMath::Cross(ZVec,axis,orthoAxis);
        double normOrthoAxis( vtkMath::Norm(orthoAxis) );
        if(normOrthoAxis > 1e-5)
        {
          //std::cerr << "ortho : " << normOrthoAxis << " X = " << orthoAxis[0] << " Y = " << orthoAxis[1] << " Z = " << orthoAxis[2] << std::endl;
          orthoAxis[0] *= normOrthoAxis; orthoAxis[1] *= normOrthoAxis; orthoAxis[2] *= normOrthoAxis;
          const double Center[3] = {0.,0.,0.};
          vtkNew<vtkDoubleArray> newArray;
          newArray->SetNumberOfComponents(3); newArray->SetNumberOfTuples(nbPoints);
          vtkDoubleArray *oldPts( vtkDoubleArray::SafeDownCast( outputMesh->GetPoints()->GetData() ) );
          double angle(asin(normOrthoAxis));
          if( vtkMath::Dot(ZVec,axis) < 0. )
            angle = M_PI - angle;
          Rotate3DAlg<double>(Center,orthoAxis,angle,nbPoints,oldPts->GetPointer(0),newArray->GetPointer(0));
          outputMesh->GetPoints()->SetData(newArray);
        }
      }
      // Translation
      {
        vtkDoubleArray *coords(vtkDoubleArray::SafeDownCast( outputMesh->GetPoints()->GetData()));
        vtkIdType nbPts(coords->GetNumberOfTuples());
        double *coordsPtr(coords->GetPointer(0));
        for(vtkIdType i = 0 ; i < nbPts ; ++i)
        { coordsPtr[3*i] += center[0]; coordsPtr[3*i+1] += center[1]; coordsPtr[3*i+2] += center[2]; }
      }
      //
      vtkNew<vtkResampleWithDataSet> probeFilter;
      probeFilter->SetInputData(outputMesh);
      probeFilter->SetSourceData(usgIn);
      //
      vtkNew<vtkPointDataToCellData> pd2cd;
      pd2cd->SetInputConnection(probeFilter->GetOutputPort());
      pd2cd->Update();
      //
      vtkNew<vtkMeshQuality> mq;
      mq->SetInputData(outputMesh);
      mq->SetTriangleQualityMeasureToArea();
      mq->SetQuadQualityMeasureToArea();
      mq->Update();
      double *area( ( vtkDoubleArray::SafeDownCast(mq->GetOutput()->GetCellData()->GetArray("Quality")) )->GetPointer(0) );
      //
      vtkInformation *outInfo(outputVector->GetInformationObject(0));
      vtkUnstructuredGrid *output(vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
      output->ShallowCopy(outputMesh);
      //
      vtkIdType nbOfCells(output->GetNumberOfCells());
      vtkFieldData* dsa(pd2cd->GetOutput()->GetCellData());
      int nbOfArrays(dsa->GetNumberOfArrays());
      for(int i = 0 ; i < nbOfArrays ; ++i )
      {
        vtkDoubleArray *arr( vtkDoubleArray::SafeDownCast(dsa->GetArray(i)) );
        if( arr && arr->GetNumberOfComponents() == 3 )
        {
          vtkNew<vtkDoubleArray> arr2;
          arr2->ShallowCopy(arr);
          output->GetCellData()->AddArray(arr2);
          double flux(ComputeFlux(nbOfCells,area,arr->GetPointer(0),axis));
          std::ostringstream oss; oss << dsa->GetArrayName(i) << "_flux";
          vtkNew<vtkDoubleArray> arrFlux;
          arrFlux->SetName(oss.str().c_str());
          arrFlux->SetNumberOfComponents(1); arrFlux->SetNumberOfTuples(nbOfCells);
          std::for_each(arrFlux->GetPointer(0),arrFlux->GetPointer(nbOfCells),[flux](double& elt) { elt = flux; });
          output->GetCellData()->AddArray(arrFlux);
        }
      }
    }
  catch(MZCException& e)
    {
      vtkErrorMacro("Exception has been thrown in vtkElectromagnetismFluxDisc::RequestInformation : " << e.what());
      return 0;
    }
  return 1;
}

void vtkElectromagnetismFluxDisc::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkElectromagnetismFluxDisc::SetCutFunction(vtkImplicitFunction* func)
{
  vtkCylinder *cyl(vtkCylinder::SafeDownCast(func));
  if(cyl)
    {
      _cyl=cyl;
      this->Modified();
    }
}

vtkMTimeType vtkElectromagnetismFluxDisc::GetMTime()
{
  vtkMTimeType maxMTime = this->Superclass::GetMTime(); // My MTime
  if(_cyl)
    {
      maxMTime=std::max(maxMTime,_cyl->GetMTime());
    }
  return maxMTime;
}
