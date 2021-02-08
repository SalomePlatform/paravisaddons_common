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

#include "vtkRateOfFlowThroughSection.h"
#include "vtkExplodePolyLine.h"
#include <vtkAdjacentVertexIterator.h>
#include <vtkAlgorithmOutput.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCharArray.h>
#include <vtkDataArraySelection.h>
#include <vtkDataObjectTreeIterator.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkDoubleArray.h>
#include <vtkInEdgeIterator.h>
#include <vtkInformation.h>
#include <vtkInformationDataObjectKey.h>
#include <vtkInformationStringKey.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkResampleWithDataSet.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkTable.h>
#include <vtkTimeStamp.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVariantArray.h>
#include <vtkWarpScalar.h>

#include "VTKToMEDMem.h"

#include <map>
#include <deque>
#include <sstream>

vtkStandardNewMacro(vtkRateOfFlowThroughSection);

///////////////////

static vtkDataSet *SplitSingleMultiBloc(vtkDataObject *ds)
{
  if(!ds)
    throw INTERP_KERNEL::Exception("vtkSedimentDeposit  SplitSingleMultiBloc : nullptr !");
  vtkMultiBlockDataSet *ds0(vtkMultiBlockDataSet::SafeDownCast(ds));
  if(!ds0)
  {
    vtkDataSet *ds00(vtkDataSet::SafeDownCast(ds));
    if(!ds00)
      throw INTERP_KERNEL::Exception("vtkSedimentDeposit  SplitSingleMultiBloc : neither a vtkMultiBlockDataSet nor a vtkDataSet !");
    return ds00;
  }
  if(ds0->GetNumberOfBlocks() != 1)
  {
    std::ostringstream oss; oss << "vtkSedimentDeposit  SplitSingleMultiBloc : presence of multiblock dataset with not exactly one dataset in it ! (" << ds0->GetNumberOfBlocks() << ") !";
    throw INTERP_KERNEL::Exception(oss.str());
  }
  vtkDataObject *ds1(ds0->GetBlock(0));
  vtkDataSet *ds1c(vtkDataSet::SafeDownCast(ds1));
  if(!ds1c)
    throw INTERP_KERNEL::Exception("vtkSedimentDeposit  SplitSingleMultiBloc : nullptr inside single multiblock element !");
  return ds1c;
}

static void ExtractInfo(vtkInformationVector *inputVector, vtkUnstructuredGrid *&usgIn)
{
  vtkInformation *inputInfo(inputVector->GetInformationObject(0));
  vtkDataSet *input(0);
  vtkDataSet *input0(vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  vtkMultiBlockDataSet *input1(vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  if (input0)
    input = input0;
  else
  {
    if (!input1)
      throw INTERP_KERNEL::Exception("Input dataSet must be a DataSet or single elt multi block dataset expected !");
    if (input1->GetNumberOfBlocks() != 1)
      throw INTERP_KERNEL::Exception("Input dataSet is a multiblock dataset with not exactly one block ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    vtkDataObject *input2(input1->GetBlock(0));
    if (!input2)
      throw INTERP_KERNEL::Exception("Input dataSet is a multiblock dataset with exactly one block but this single element is NULL !");
    vtkDataSet *input2c(vtkDataSet::SafeDownCast(input2));
    if (!input2c)
      throw INTERP_KERNEL::Exception("Input dataSet is a multiblock dataset with exactly one block but this single element is not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    input = input2c;
  }
  if (!input)
    throw INTERP_KERNEL::Exception("Input data set is NULL !");
  usgIn = vtkUnstructuredGrid::SafeDownCast(input);
  if (!usgIn)
    throw INTERP_KERNEL::Exception("Input data set is not an unstructured mesh ! This filter works only on unstructured meshes !");
}

////////////////////

void vtkRateOfFlowThroughSection::vtkInternal::fillTable(vtkTable *table) const
{
  {
    vtkNew<vtkDoubleArray> timeArr;
    timeArr->SetName("Time");
    timeArr->SetNumberOfTuples(_data.size());
    double *pt(timeArr->GetPointer(0));
    {
      std::size_t tmp(0);
      std::for_each(pt, pt + _data.size(), [this, &tmp](double &val) { val = this->_data[tmp++].first; });
    }
    table->AddColumn(timeArr);
  }
  {
    vtkNew<vtkDoubleArray> timeArr;
    timeArr->SetName("Rate of flow");
    timeArr->SetNumberOfTuples(_data.size());
    double *pt(timeArr->GetPointer(0));
    {
      std::size_t tmp(0);
      std::for_each(pt, pt + _data.size(), [this, &tmp](double &val) { val = this->_data[tmp++].second; });
    }
    table->AddColumn(timeArr);
  }
}

void vtkRateOfFlowThroughSection::vtkInternal::analyzeInputDataSets(vtkUnstructuredGrid *ds1, vtkDataSet *ds2)
{
  _recomputationOfMatrixNeeded = false;
  if (_mt1 != ds1->GetMeshMTime())
  {
    _mt1 = ds1->GetMeshMTime();
    _recomputationOfMatrixNeeded = true;
  }
  vtkUnstructuredGrid *ds2_0(vtkUnstructuredGrid::SafeDownCast(ds2));
  vtkPolyData *ds2_1(vtkPolyData::SafeDownCast(ds2));
  if (!ds2_0 && !ds2_1)
    throw INTERP_KERNEL::Exception("analyzeInputDataSets : unexpected source !");
  if (ds2_0)
    if (_mt2 != ds2_0->GetMeshMTime())
    {
      _mt2 = ds2_0->GetMeshMTime();
      _recomputationOfMatrixNeeded = true;
    }
  if (ds2_1)
    if (_mt2 != ds2_1->GetMeshMTime())
    {
      _mt2 = ds2_1->GetMeshMTime();
      _recomputationOfMatrixNeeded = true;
    }
}

////////////////////

vtkRateOfFlowThroughSection::vtkRateOfFlowThroughSection() : NumberOfTimeSteps(0), CurrentTimeIndex(0), IsExecuting(false), Internal(nullptr)
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

vtkRateOfFlowThroughSection::~vtkRateOfFlowThroughSection()
{
}

int vtkRateOfFlowThroughSection::RequestUpdateExtent(vtkInformation *, vtkInformationVector **inputVector, vtkInformationVector *vtkNotUsed(outputVector))
{
  // vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkInformation *inInfo1 = inputVector[0]->GetInformationObject(0);

  // get the requested update extent
  double *inTimes = inInfo1->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  if (inTimes)
  {
    double timeReq = inTimes[this->CurrentTimeIndex];
    inInfo1->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), timeReq);
  }

  return 1;
}

int vtkRateOfFlowThroughSection::RequestInformation(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkRateOfFlowThroughSection::RequestInformation ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid *usgIn(0);
    ExtractInfo(inputVector[0], usgIn);
    vtkInformation *info(outputVector->GetInformationObject(0));
    vtkInformation *inInfo(inputVector[0]->GetInformationObject(0));
    if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
    {
      this->NumberOfTimeSteps = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    }
    else
    {
      this->NumberOfTimeSteps = 0;
    }
    // The output of this filter does not contain a specific time, rather
    // it contains a collection of time steps. Also, this filter does not
    // respond to time requests. Therefore, we remove all time information
    // from the output.
    vtkInformation *outInfo(outputVector->GetInformationObject(0));
    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
    {
      outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    }
    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_RANGE()))
    {
      outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());
    }
  }
  catch (INTERP_KERNEL::Exception &e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkRateOfFlowThroughSection::RequestInformation : " << e.what() << std::endl;
    if (this->HasObserver("ErrorEvent"))
      this->InvokeEvent("ErrorEvent", const_cast<char *>(oss.str().c_str()));
    else
      vtkOutputWindowDisplayErrorText(const_cast<char *>(oss.str().c_str()));
    vtkObject::BreakOnError();
    return 0;
  }
  return 1;
}

static MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> ToMedcoupling(MEDCoupling::MCAuto<MEDCoupling::MEDFileData> &mfd, vtkDataSet *usgIn)
{
  WriteMEDFileFromVTKDataSet(mfd, usgIn, {}, 0., 0);
  MEDCoupling::MEDFileMeshes *ms(mfd->getMeshes());
  if (ms->getNumberOfMeshes() != 1)
    throw INTERP_KERNEL::Exception("Unexpected number of meshes !");
  MEDCoupling::MEDFileMesh *mm(ms->getMeshAtPos(0));
  MEDCoupling::MEDFileUMesh *mmu(dynamic_cast<MEDCoupling::MEDFileUMesh *>(mm));
  if (!mmu)
    throw INTERP_KERNEL::Exception("Expecting unstructured one !");
  return mmu->getMeshAtLevel(0);
}

static void MyAssert(bool status, const std::string &message)
{
  if (!status)
    throw INTERP_KERNEL::Exception(message);
}

bool IsNameIn(const std::string& name, const std::vector<std::string>& namesPossible)
{
  for(auto np : namesPossible)
    {
      std::size_t pos( name.find(np) );
      if(pos==std::string::npos)
        continue;
      std::string nameCpy(name);
      std::string tmp(nameCpy.replace(pos,np.length(),std::string()));
      if( tmp.find_first_not_of(" \t") == std::string::npos )
        return true;
    }
  return false;
}

vtkDataArray *FindArrayHavingNameIn(vtkPointData *pd, const std::vector<std::string>& namesPossible, std::function<bool(vtkDataArray *)> func)
{
  vtkDataArray *ret(nullptr);
  for(auto i = 0; i < pd->GetNumberOfArrays() ; ++i )
    {
      vtkDataArray *arr(pd->GetArray(i));
      std::string name(arr->GetName());
      if(IsNameIn(name,namesPossible))
        {
          if( func(arr) )
            {
              ret = arr;
              break;
            }
        }
    }
  return ret;
}

int vtkRateOfFlowThroughSection::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkRateOfFlowThroughSection::RequestData        ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid *usgIn(nullptr);
    ExtractInfo(inputVector[0], usgIn);
    // is this the first request
    if (!this->IsExecuting)
    {
      request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
      this->IsExecuting = true;
      delete this->Internal;
      this->Internal = new vtkInternal;
    }
    //
    vtkInformation *sourceInfo(inputVector[1]->GetInformationObject(0));
    vtkDataObject *source(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkDataSet *source1(SplitSingleMultiBloc(source));
    //
    vtkNew<vtkExplodePolyLine> epl;
    epl->SetInputData(source1);
    epl->Update();
    vtkDataSet *source2(epl->GetOutput());
    //
    this->Internal->analyzeInputDataSets(usgIn, source1);
    ///////////////////////
    //////////////////////
    /////////////////////
    //
    std::vector<std::map<int, double>> &matrix(this->Internal->getMatrix());
    if (this->Internal->computationNeeded())
    {
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> m, sec;
      MEDCoupling::MCAuto<MEDCoupling::MEDFileData> mfd(MEDCoupling::MEDFileData::New());
      m = ToMedcoupling(mfd, usgIn);
      {
        MEDCoupling::MCAuto<MEDCoupling::MEDFileData> mfdSec(MEDCoupling::MEDFileData::New());
        sec = ToMedcoupling(mfd, source2);
      }
      {
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> arr(m->getCoords()->keepSelectedComponents({2}));
        MyAssert(arr->isUniform(0, 1e-12), "Expected coords array equal to 0 for Z axis.");
      }
      m->changeSpaceDimension(2, 0.);
      {
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> arr(sec->getCoords()->keepSelectedComponents({2}));
        MyAssert(arr->isUniform(0, 1e-12), "Expected coords array equal to 0 for Z axis.");
      }
      sec->changeSpaceDimension(2, 0.);
      // sec peut etre completement merdique avec des pts dupliques alors on filtre
      {
        bool tmp;
        mcIdType tmpp;
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> tmppp(sec->mergeNodes(1e-12, tmp, tmpp));
      }
      sec->zipCoords();
      sec->removeDegenerated1DCells();
      //
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> line_inter;
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> cellid_in_2d, cellid_in1d;
      {
        MEDCoupling::MEDCouplingUMesh *tmp(nullptr), *tmp2(nullptr);
        MEDCoupling::DataArrayIdType *tmp3(nullptr), *tmp4(nullptr);
        MEDCoupling::MEDCouplingUMesh::Intersect2DMeshWith1DLine(m, sec, 1e-12, tmp, tmp2, tmp3, tmp4);
        tmp->decrRef();
        line_inter = tmp2;
        cellid_in_2d = tmp3;
        cellid_in1d = tmp4;
      }
      line_inter->zipCoords();
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> TwoDcells(MEDCoupling::DataArrayIdType::New());
      TwoDcells->alloc(line_inter->getNumberOfCells(), 1);
      {
        auto t(cellid_in1d->begin());
        auto TwoDcellsPtr(TwoDcells->getPointer());
        for (std::size_t i = 0; i < cellid_in1d->getNumberOfTuples(); i++, t += 2, TwoDcellsPtr++)
        {
          int zeValue(-1);
          std::for_each(t, t + 2, [&zeValue](const int &v) { if(v!=-1 && zeValue==-1) zeValue=v; });
          *TwoDcellsPtr = zeValue;
        }
      }
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> notFreeStyle1DCells(TwoDcells->findIdsNotEqual(-1));
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> n2oCells(TwoDcells->selectByTupleId(notFreeStyle1DCells->begin(), notFreeStyle1DCells->end()));
      TwoDcells = cellid_in_2d->selectByTupleId(n2oCells->begin(), n2oCells->end());
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> effective_line1d(line_inter->buildPartOfMySelf(notFreeStyle1DCells->begin(), notFreeStyle1DCells->end()));
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> effective_2d_cells(m->buildPartOfMySelf(TwoDcells->begin(), TwoDcells->end()));
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> o2n(effective_2d_cells->zipCoordsTraducer());
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> n2o(o2n->invertArrayO2N2N2O(effective_2d_cells->getNumberOfNodes()));
      MEDCoupling::MCAuto<MEDCoupling::MEDCoupling1SGTUMesh> effective_line1d_2(MEDCoupling::MEDCoupling1SGTUMesh::New(effective_line1d));     // change format of umesh to ease alg
      MEDCoupling::MCAuto<MEDCoupling::MEDCoupling1SGTUMesh> effective_2d_cells_2(MEDCoupling::MEDCoupling1SGTUMesh::New(effective_2d_cells)); // change format of umesh to ease alg
      MyAssert(effective_2d_cells_2->getCellModelEnum() == INTERP_KERNEL::NORM_TRI3, "Only TRI3 are expected as input");

      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> conn1d(effective_line1d_2->getNodalConnectivity()->deepCopy());
      conn1d->rearrange(2);
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> conn2d(effective_2d_cells_2->getNodalConnectivity()->deepCopy());
      conn2d->rearrange(3);
      const MEDCoupling::DataArrayDouble *coo1d(effective_line1d->getCoords()), *coo2d(effective_2d_cells->getCoords());
      MyAssert(conn2d->getNumberOfTuples() == conn1d->getNumberOfTuples(), "Internal error !");

      matrix.resize(effective_line1d->getNumberOfCells());
      {
        const mcIdType *t1(conn1d->begin()), *t2(conn2d->begin());
        const double *coo1dPtr(coo1d->begin()), *coo2dPtr(coo2d->begin());
        double seg2[4], tri3[6];
        for (std::size_t i = 0; i < effective_line1d->getNumberOfCells(); i++, t1 += 2, t2 += 3)
        {
          double baryInfo[3];
          double length;
          seg2[0] = coo1dPtr[2 * t1[0]];
          seg2[1] = coo1dPtr[2 * t1[0] + 1];
          seg2[2] = coo1dPtr[2 * t1[1]];
          seg2[3] = coo1dPtr[2 * t1[1] + 1];
          tri3[0] = coo2dPtr[2 * t2[0]];
          tri3[1] = coo2dPtr[2 * t2[0] + 1];
          tri3[2] = coo2dPtr[2 * t2[1]];
          tri3[3] = coo2dPtr[2 * t2[1] + 1];
          tri3[4] = coo2dPtr[2 * t2[2]];
          tri3[5] = coo2dPtr[2 * t2[2] + 1];
          MEDCoupling::DataArrayDouble::ComputeIntegralOfSeg2IntoTri3(seg2, tri3, baryInfo, length);
          std::map<int, double> &row(matrix[i]);
          row[n2o->getIJ(t2[0], 0)] = baryInfo[0];
          row[n2o->getIJ(t2[1], 0)] = baryInfo[1];
          row[n2o->getIJ(t2[2], 0)] = baryInfo[2];
        }
      }
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> orthoField;
      const MEDCoupling::DataArrayDouble *ortho(nullptr);
      {
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> tmp(effective_line1d->buildUnstructured());
        orthoField = tmp->buildOrthogonalField();
        this->Internal->setOrtho(orthoField->getArray());
      }
      {
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> measure(effective_line1d->getMeasureField(true));
        this->Internal->setMeasure(measure->getArray());
      }
    }
    const MEDCoupling::DataArrayDouble *ortho(this->Internal->getOrtho()), *measure_arr(this->Internal->getMeasure());
    ///////////////////////
    //////////////////////
    /////////////////////
    constexpr char SEARCHED_FIELD_HAUTEUR[]="HAUTEUR D'EAU";
    constexpr char SEARCHED_FIELD_HAUTEUR2[]="WATER DEPTH";
    constexpr char SEARCHED_FIELD_SPEED[]="VITESSE *";
    constexpr char SEARCHED_FIELD_SPEED2[]="VELOCITY *";

    vtkPointData *pd(usgIn->GetPointData());
    vtkDataArray *h_water_tmp(FindArrayHavingNameIn(pd,{ SEARCHED_FIELD_HAUTEUR, SEARCHED_FIELD_HAUTEUR2 },[](vtkDataArray *arr) { return arr->GetNumberOfComponents()==1; }));
    vtkDataArray *speed_tmp(FindArrayHavingNameIn(pd,{ "VITESSE *", "VELOCITY *" },[](vtkDataArray *arr) { return arr->GetNumberOfComponents()>1; }));
    vtkDoubleArray *h_water(vtkDoubleArray::SafeDownCast(h_water_tmp)), *speed(vtkDoubleArray::SafeDownCast(speed_tmp));
    std::ostringstream oss;
    oss << "Expecting presence of float32 following fields : \"" << SEARCHED_FIELD_HAUTEUR << "\" or \"" << SEARCHED_FIELD_HAUTEUR2 << "\"and \"" << SEARCHED_FIELD_SPEED << "\" or \"" << SEARCHED_FIELD_SPEED2 << " \"";
    MyAssert(h_water && speed, oss.str());
    MEDCoupling::MCAuto<MEDCoupling::DataArrayFloat> speed2(MEDCoupling::DataArrayFloat::New());
    speed2->alloc(speed->GetNumberOfTuples(), speed->GetNumberOfComponents());
    std::copy(speed->GetPointer(0), speed->GetPointer(0) + speed2->getNbOfElems(), speed2->rwBegin());
    MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> speed_arr(speed2->convertToDblArr());
    MEDCoupling::MCAuto<MEDCoupling::DataArrayFloat> h_water_arrf(MEDCoupling::DataArrayFloat::New());
    h_water_arrf->alloc(h_water->GetNumberOfTuples(), 1);
    std::copy(h_water->GetPointer(0), h_water->GetPointer(0) + h_water_arrf->getNbOfElems(), h_water_arrf->rwBegin());
    /*MEDCoupling::MEDFileFloatField1TS *speed(nullptr),*h_water(nullptr);
      {
        MEDCoupling::MEDFileFields *fields(mfd->getFields());
        MyAssert(fields,"No arrays found in the input dataset !");
        MEDCoupling::MEDFileAnyTypeFieldMultiTS *speed_mts(fields->getFieldWithName(SEARCHED_FIELD_SPEED));
        MEDCoupling::MEDFileAnyTypeFieldMultiTS *h_water_mts(fields->getFieldWithName(SEARCHED_FIELD_HAUTEUR));
        std::ostringstream oss; oss << "Expecting single time step for following fields : \"" << SEARCHED_FIELD_HAUTEUR << "\" and \"" <<  SEARCHED_FIELD_SPEED << "\"";
        MyAssert(speed_mts->getNumberOfTS()==1 && h_water_mts->getNumberOfTS()==1,oss.str());
        MEDCoupling::MEDFileFloatFieldMultiTS *speed_mts_2(dynamic_cast<MEDCoupling::MEDFileFloatFieldMultiTS *>(speed_mts));
        MEDCoupling::MEDFileFloatFieldMultiTS *h_water_mts_2(dynamic_cast<MEDCoupling::MEDFileFloatFieldMultiTS *>(h_water_mts));
        std::ostringstream oss2; oss2 << "Expecting presence of float32 following fields : \"" << SEARCHED_FIELD_HAUTEUR << "\" and \"" <<  SEARCHED_FIELD_SPEED << "\"";
        MyAssert(!speed_mts_2 || !h_water_mts_2,oss2.str());
        speed=speed_mts_2->getTimeStepAtPos(0); h_water=h_water_mts_2->getTimeStepAtPos(0);
      }
      MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> speed_arr(speed->getUndergroundDataArray()->convertToDblArr());*/
    {
      MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> arr(speed_arr->keepSelectedComponents({2}));
      MyAssert(arr->isUniform(0, 1e-12), "Expected speed array equal to 0 for Z axis.");
    }
    speed_arr = speed_arr->keepSelectedComponents({0, 1});
    MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> h_water_arr(h_water_arrf->convertToDblArr());
    const double *speed_ptr(speed_arr->begin()), *ortho_ptr(ortho->begin());
    MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> h_out(MEDCoupling::DataArrayDouble::New());
    h_out->alloc(matrix.size(), 1);
    h_out->fillWithZero();
    MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> v_out(MEDCoupling::DataArrayDouble::New());
    v_out->alloc(matrix.size(), 1);
    v_out->fillWithZero();
    const double *orthoPtr(ortho->begin());
    for (std::size_t i = 0; i < matrix.size(); i++, speed_ptr += 2, ortho_ptr += 2)
    {
      const std::map<int, double> &row(matrix[i]);
      double h_out_value(0.), speed[2] = {0., 0.};
      for (auto it : row)
      {
        h_out_value += it.second * h_water_arr->getIJ(it.first, 0);
        speed[0] += it.second * speed_arr->getIJ(it.first, 0);
        speed[1] += it.second * speed_arr->getIJ(it.first, 1);
      }
      h_out->setIJ(i, 0, h_out_value);
      v_out->setIJ(i, 0, speed[0]*orthoPtr[2*i+0] +speed[1]*orthoPtr[2*i+1] );
    }
    double zeValue(0.);
    {
      MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp1(MEDCoupling::DataArrayDouble::Multiply(h_out, v_out));
      tmp1 = MEDCoupling::DataArrayDouble::Multiply(tmp1, measure_arr);
      zeValue = std::abs(tmp1->accumulate((std::size_t)0));
    }
    double timeStep;
    {
      vtkInformation *inInfo(inputVector[0]->GetInformationObject(0));
      vtkDataObject *input(vtkDataObject::GetData(inInfo));
      timeStep = input->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());
    }
    this->Internal->pushData(timeStep, zeValue);
    this->CurrentTimeIndex++;
    if (this->CurrentTimeIndex == this->NumberOfTimeSteps)
    {
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      this->CurrentTimeIndex = 0;
      this->IsExecuting = false;
      vtkInformation *outInfo(outputVector->GetInformationObject(0));
      vtkTable *output(vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
      vtkNew<vtkTable> table;
      this->Internal->fillTable(table);
      output->ShallowCopy(table);
    }
  }
  catch (INTERP_KERNEL::Exception& e)
  {
    if (this->IsExecuting)
    {
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      this->CurrentTimeIndex = 0;
      this->IsExecuting = false;
    }
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkRateOfFlowThroughSection::RequestData : " << e.what() << std::endl;
    vtkErrorMacro(<< oss.str());
    return 0;
  }
  return 1;
}

void vtkRateOfFlowThroughSection::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkRateOfFlowThroughSection::SetSourceData(vtkDataObject *input)
{
  this->SetInputData(1, input);
}

void vtkRateOfFlowThroughSection::SetSourceConnection(vtkAlgorithmOutput *algOutput)
{
  this->SetInputConnection(1, algOutput);
}

int vtkRateOfFlowThroughSection::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
  return 1;
}
