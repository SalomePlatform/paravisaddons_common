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
// Author : Anthony Geay

#include "vtkElectromagnetismRotation.h"

#include "MEDCouplingMemArray.hxx"

#include "ElectromagnetismRotationHelper.h"
#include "VTKMEDTraits.hxx"

#include "vtkAdjacentVertexIterator.h"
#include "vtkAOSDataArrayTemplate.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#ifdef WIN32
#include "vtkLongLongArray.h"
#endif
#include "vtkCellData.h"
#include "vtkPointData.h"

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
#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkDemandDrivenPipeline.h"
#include "vtkDataObjectTreeIterator.h"
#include "vtkThreshold.h"
#include "vtkMultiBlockDataGroupFilter.h"
#include "vtkMergeBlocks.h"
#include "vtkInformationDataObjectMetaDataKey.h"
#include "vtkTransform.h"
#include "vtkFunctionParser.h"

#include <map>
#include <deque>
#include <cstring>
#include <memory>

const char ZE_SEP[]="@@][@@";

const char TS_STR[]="TS";

const char COM_SUP_STR[]="ComSup";

const char FAMILY_ID_CELL_NAME[]="FamilyIdCell";

const char NUM_ID_CELL_NAME[]="NumIdCell";

const char FAMILY_ID_NODE_NAME[]="FamilyIdNode";

const char NUM_ID_NODE_NAME[]="NumIdNode";

const char GLOBAL_NODE_ID_NAME[]="GlobalNodeIds";

vtkStandardNewMacro(vtkElectromagnetismRotation);

static vtkInformationDataObjectMetaDataKey* GetMEDReaderMetaDataIfAny()
{
  static const char ZE_KEY[] = "vtkFileSeriesGroupReader::META_DATA";
  MEDCoupling::GlobalDict* gd(MEDCoupling::GlobalDict::GetInstance());
  if (!gd->hasKey(ZE_KEY))
    return 0;
  std::string ptSt(gd->value(ZE_KEY));
  void* pt(0);
  std::istringstream iss(ptSt);
  iss >> pt;
  return reinterpret_cast<vtkInformationDataObjectMetaDataKey*>(pt);
}

bool IsInformationOK(vtkInformation* info)
{
  vtkInformationDataObjectMetaDataKey* key(GetMEDReaderMetaDataIfAny());
  if (!key)
    return false;
  // Check the information contain meta data key
  if (!info->Has(key))
    return false;
  // Recover Meta Data
  vtkMutableDirectedGraph* sil(vtkMutableDirectedGraph::SafeDownCast(info->Get(key)));
  if (!sil)
    return false;
  int idNames(0);
  vtkAbstractArray* verticesNames(sil->GetVertexData()->GetAbstractArray("Names", idNames));
  vtkStringArray* verticesNames2(vtkStringArray::SafeDownCast(verticesNames));
  if (!verticesNames2)
    return false;
  for (int i = 0; i < verticesNames2->GetNumberOfValues(); i++)
  {
    vtkStdString& st(verticesNames2->GetValue(i));
    if (st == "MeshesFamsGrps")
      return true;
  }
  return false;
}

class vtkElectromagnetismRotation::vtkElectromagnetismRotationInternal : public ElectromagnetismRotationInternal
{
};

////////////////////

vtkElectromagnetismRotation::vtkElectromagnetismRotation():SIL(NULL),Internal(new vtkElectromagnetismRotationInternal),InsideOut(0),Axis(2),RotationRotor(1e300)
{
}

vtkElectromagnetismRotation::~vtkElectromagnetismRotation()
{
  delete this->Internal;
}

void vtkElectromagnetismRotation::SetInsideOut(int val)
{
  if(this->InsideOut!=val)
    {
      this->InsideOut=val;
      this->Modified();
    }
}

void vtkElectromagnetismRotation::SetAxis(int axis)
{
  if(this->Axis!=axis)
    {
      this->Axis=axis;
      this->Modified();
    }
}

void vtkElectromagnetismRotation::SetAngularStep(char *angStep)
{
  this->AngularStep = angStep;
  this->Modified();
}

int vtkElectromagnetismRotation::RequestInformation(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
//  vtkUnstructuredGridAlgorithm::RequestInformation(request,inputVector,outputVector);
  try
    {
//      std::cerr << "########################################## vtkElectromagnetismRotation::RequestInformation ##########################################" << std::endl;
//      request->Print(cout);
      vtkInformation *outInfo(outputVector->GetInformationObject(0));
      vtkInformation *inputInfo(inputVector[0]->GetInformationObject(0));
      if(!ElectromagnetismRotationInternal::IndependantIsInformationOK(GetMEDReaderMetaDataIfAny(),inputInfo))
        {
        vtkErrorMacro("No SIL Data available ! The source of this filter must be MEDReader !");
        return 0;
        }

      this->SetSIL(vtkMutableDirectedGraph::SafeDownCast(inputInfo->Get(GetMEDReaderMetaDataIfAny())));
      this->Internal->loadFrom(this->SIL);
      //this->Internal->printMySelf(std::cerr);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      std::cerr << "Exception has been thrown in vtkElectromagnetismRotation::RequestInformation : " << e.what() << std::endl;
      return 0;
    }
  return 1;
}

/*!
 * Do not use vtkCxxSetObjectMacro macro because input mdg comes from an already managed in the pipeline just a ref on it.
 */
void vtkElectromagnetismRotation::SetSIL(vtkMutableDirectedGraph *mdg)
{
  if(this->SIL==mdg)
    return ;
  this->SIL=mdg;
}

template<class CellPointExtractor>
vtkDataSet *FilterFamilies(vtkThreshold *thres,
                           vtkDataSet *input, const std::set<int>& idsToKeep, bool insideOut, const char *arrNameOfFamilyField,
                           const char *associationForThreshold, bool& catchAll, bool& catchSmth)
{
  const int VTK_DATA_ARRAY_DELETE=vtkAOSDataArrayTemplate<double>::VTK_DATA_ARRAY_DELETE;
  const char ZE_SELECTION_ARR_NAME[]="@@ZeSelection@@";
  vtkDataSet *output(input->NewInstance());
  output->ShallowCopy(input);
  thres->SetInputData(output);
  vtkDataSetAttributes *dscIn(input->GetCellData()),*dscIn2(input->GetPointData());
  vtkDataSetAttributes *dscOut(output->GetCellData()),*dscOut2(output->GetPointData());
  //
  double vMin(insideOut==0?1.:0.),vMax(insideOut==0?2.:1.);
  thres->SetLowerThreshold(vMin);
  thres->SetUpperThreshold(vMax);
  // OK for the output
  //
  CellPointExtractor cpe2(input);
  vtkDataArray *da(cpe2.Get()->GetScalars(arrNameOfFamilyField));
  if(!da)
    return 0;
  std::string daName(da->GetName());
  typedef MEDFileVTKTraits<mcIdType>::VtkType vtkMCIdTypeArray;
  vtkMCIdTypeArray *dai(vtkMCIdTypeArray::SafeDownCast(da));
  if(daName!=arrNameOfFamilyField || !dai)
    return 0;
  //
  int nbOfTuples(dai->GetNumberOfTuples());
  vtkCharArray *zeSelection(vtkCharArray::New());
  zeSelection->SetName(ZE_SELECTION_ARR_NAME);
  zeSelection->SetNumberOfComponents(1);
  char *pt(new char[nbOfTuples]);
  zeSelection->SetArray(pt,nbOfTuples,0,VTK_DATA_ARRAY_DELETE);
  const mcIdType *inPtr(dai->GetPointer(0));
  std::fill(pt,pt+nbOfTuples,0);
  catchAll=true; catchSmth=false;
  std::vector<bool> pt2(nbOfTuples,false);
  for(std::set<int>::const_iterator it=idsToKeep.begin();it!=idsToKeep.end();it++)
    {
      bool catchFid(false);
      for(int i=0;i<nbOfTuples;i++)
        if(inPtr[i]==*it)
          { pt2[i]=true; catchFid=true; }
      if(!catchFid)
        catchAll=false;
      else
        catchSmth=true;
    }
  for(int ii=0;ii<nbOfTuples;ii++)
    if(pt2[ii])
      pt[ii]=2;
  CellPointExtractor cpe3(output);
  int idx(cpe3.Get()->AddArray(zeSelection));
  cpe3.Get()->SetActiveAttribute(idx,vtkDataSetAttributes::SCALARS);
  cpe3.Get()->CopyScalarsOff();
  zeSelection->Delete();
  //
  thres->SetInputArrayToProcess(idx,0,0,associationForThreshold,ZE_SELECTION_ARR_NAME);
  thres->Update();
  vtkUnstructuredGrid *zeComputedOutput(thres->GetOutput());
  CellPointExtractor cpe(zeComputedOutput);
  cpe.Get()->RemoveArray(idx);
  output->Delete();
  zeComputedOutput->Register(0);
  return zeComputedOutput;
}

class CellExtractor
{
public:
  CellExtractor(vtkDataSet *ds):_ds(ds) { }
  vtkDataSetAttributes *Get() { return _ds->GetCellData(); }
private:
  vtkDataSet *_ds;
};

class PointExtractor
{
public:
  PointExtractor(vtkDataSet *ds):_ds(ds) { }
  vtkDataSetAttributes *Get() { return _ds->GetPointData(); }
private:
  vtkDataSet *_ds;
};

int vtkElectromagnetismRotation::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  try
    {
      // std::cerr << "########################################## vtkElectromagnetismRotation::RequestData        ##########################################" << std::endl;
      // request->Print(cout);
      vtkInformation* inputInfo=inputVector[0]->GetInformationObject(0);
      vtkMultiBlockDataSet *inputMB(vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
      if(inputMB->GetNumberOfBlocks()!=1)
        {
          vtkErrorMacro(<< "vtkElectromagnetismRotation::RequestData : input has not the right number of parts ! Expected 1 !" ) ;
          return 0;
        }
      vtkDataSet *input(vtkDataSet::SafeDownCast(inputMB->GetBlock(0)));
      vtkInformation *info(input->GetInformation());
      vtkInformation *outInfo(outputVector->GetInformationObject(0));
      vtkMultiBlockDataSet *output(vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
      output->SetNumberOfBlocks(2);
      std::set<int> idsToKeep(this->Internal->getIdsToKeep());
      this->Internal->clearSelection();
      // first shrink the input
      bool catchAll,catchSmth;
      vtkNew<vtkThreshold> thres1,thres2;
      vtkSmartPointer<vtkDataSet> rotor(FilterFamilies<CellExtractor>(thres1,input,idsToKeep,0,FAMILY_ID_CELL_NAME,"vtkDataObject::FIELD_ASSOCIATION_CELLS",catchAll,catchSmth));
      vtkSmartPointer<vtkDataSet> stator(FilterFamilies<CellExtractor>(thres2,input,idsToKeep,1,FAMILY_ID_CELL_NAME,"vtkDataObject::FIELD_ASSOCIATION_CELLS",catchAll,catchSmth));
      //
      double reqTS(0.);
      int nbOfSteps(-1);
      std::unique_ptr<double[]> timeSteps;
       
      if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
        reqTS=outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
      if(outInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
      {
        nbOfSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        timeSteps.reset(new double[ nbOfSteps ]);
        outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),timeSteps.get());
        //std::cerr << "nb : " << nbOfSteps << std::endl;
        //std::for_each(timeSteps.get(),timeSteps.get()+nbOfSteps,[](double v) { std::cerr << v << std::endl; });
      }
      if(nbOfSteps<2 || !timeSteps.get())
      {
        vtkErrorMacro(<< "vtkElectromagnetismRotation::RequestData : A temporal dataset is expected ! Here < 2 time steps found !" ) ;
        return 0;
      }
      // Calcul de l angle effectif de rotation
      double minTime(timeSteps[0]),maxTime(timeSteps[nbOfSteps-1]);
      if(minTime == maxTime)
      {
        vtkErrorMacro(<< "vtkElectromagnetismRotation::RequestData : minTime == maxTime !" ) ;
        return 0;
      }
      double angularStepD = 0;
      {
        vtkNew<vtkFunctionParser> fp;
        fp->SetFunction(this->AngularStep.c_str());
        angularStepD = fp->GetScalarResult();
      }
      double effectiveTime(reqTS); effectiveTime = std::max(effectiveTime,minTime); effectiveTime = std::min(effectiveTime,maxTime);
      this->RotationRotor = (angularStepD*((effectiveTime-minTime)/(maxTime-minTime)));
      //std::cout << "*** " << effectiveTime << " " << minTime << " " << maxTime << " " << angleDegree << std::endl << std::flush;
      //std::cout << "*** " << this->RotationRotor << std::endl << std::flush;
      //
      if(rotor)
        {
          if(catchAll)
            {
              vtkNew<vtkTransform> transformR;
              switch(this->Axis)
              {
                case 0:
                {
                  transformR->RotateX(this->RotationRotor);
                  break;
                }
                case 1:
                {
                  transformR->RotateY(this->RotationRotor);
                  break;
                }
                case 2:
                {
                  transformR->RotateZ(this->RotationRotor);
                  break;
                }
                default:
                {
                  vtkErrorMacro(<< "vtkElectromagnetismRotation::RequestData : not recognized axis !" ) ;
                  return 0;
                }
              }
              vtkNew<vtkPoints> newCoords;
              transformR->TransformPoints(vtkPointSet::SafeDownCast(rotor)->GetPoints(),newCoords);
              vtkPointSet::SafeDownCast(rotor)->SetPoints(newCoords);
              output->SetBlock(0,rotor);
              output->SetBlock(1,stator);
              return 1;
            }
          else
            return 0;
        }
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      vtkErrorMacro(<< "Exception has been thrown in vtkElectromagnetismRotation::RequestData : " << e.what());
      return 0;
    }
}

int vtkElectromagnetismRotation::GetSILUpdateStamp()
{
  return (int)this->SILTime;
}

const char* vtkElectromagnetismRotation::GetGrpStart()
{
  return ElectromagnetismRotationGrp::start();
}

const char* vtkElectromagnetismRotation::GetFamStart()
{
  return ElectromagnetismRotationFam::start();
}

void vtkElectromagnetismRotation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

int vtkElectromagnetismRotation::GetNumberOfGroupsFlagsArrays()
{
  int ret(this->Internal->getNumberOfEntries());
  //std::cerr << "vtkElectromagnetismRotation::GetNumberOfFieldsTreeArrays() -> " << ret << std::endl;
  return ret;
}

const char *vtkElectromagnetismRotation::GetGroupsFlagsArrayName(int index)
{
  const char *ret(this->Internal->getKeyOfEntry(index));
//  std::cerr << "vtkElectromagnetismRotation::GetFieldsTreeArrayName(" << index << ") -> " << ret << std::endl;
  return ret;
}

int vtkElectromagnetismRotation::GetGroupsFlagsArrayStatus(const char *name)
{
  int ret((int)this->Internal->getStatusOfEntryStr(name));
//  std::cerr << "vtkElectromagnetismRotation::GetGroupsFlagsArrayStatus(" << name << ") -> " << ret << std::endl;
  return ret;
}

void vtkElectromagnetismRotation::SetGroupsFlagsStatus(const char *name, int status)
{
  //std::cerr << "vtkElectromagnetismRotation::SetFieldsStatus(" << name << "," << status << ")" << std::endl;
  this->Internal->setStatusOfEntryStr(name,(bool)status);
  this->Modified();
  //this->Internal->printMySelf(std::cerr);
}

const char *vtkElectromagnetismRotation::GetMeshName()
{
  return this->Internal->getMeshName();
}
