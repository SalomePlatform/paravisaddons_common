// Copyright (C) 2021-2023  CEA/DEN, EDF R&D
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

#include "vtkSedimentDeposit.h"
#include "vtkExplodePolyLine.h"
#include <vtkSetGet.h>


vtkStandardNewMacro(vtkSedimentDeposit);

///////////////////

static int GetNumberOfBlocs(vtkDataObject *ds)
{
  if(!ds)
    throw INTERP_KERNEL::Exception("vtkSedimentDeposit  SplitSingleMultiBloc : nullptr !");
  vtkMultiBlockDataSet *ds0(vtkMultiBlockDataSet::SafeDownCast(ds));
  if(!ds0)
  {
    vtkDataSet *ds00(vtkDataSet::SafeDownCast(ds));
    if(!ds00)
      throw INTERP_KERNEL::Exception("vtkSedimentDeposit  SplitSingleMultiBloc : neither a vtkMultiBlockDataSet nor a vtkDataSet !");
    return 1;
  }
  return ds0->GetNumberOfBlocks();
}

static vtkDataSet *SplitSingleMultiBloc(vtkDataObject *ds, int blockId)
{
  if(!ds)
    throw INTERP_KERNEL::Exception("vtkSedimentDeposit  SplitSingleMultiBloc : nullptr !");
  vtkMultiBlockDataSet *ds0(vtkMultiBlockDataSet::SafeDownCast(ds));
  if(!ds0)
  {
    vtkDataSet *ds00(vtkDataSet::SafeDownCast(ds));
    if(!ds00)
      throw INTERP_KERNEL::Exception("vtkSedimentDeposit  SplitSingleMultiBloc : neither a vtkMultiBlockDataSet nor a vtkDataSet !");
    if(blockId != 0)
      throw INTERP_KERNEL::Exception("vtkSedimentDeposit  SplitSingleMultiBloc : 0 expected !");
    return ds00;
  }
  if( blockId >= ds0->GetNumberOfBlocks() )
  {
    std::ostringstream oss; oss << "vtkSedimentDeposit  SplitSingleMultiBloc : presence of multiblock dataset with not exactly one dataset in it ! (" << ds0->GetNumberOfBlocks() << ") !";
    throw INTERP_KERNEL::Exception(oss.str());
  }
  vtkDataObject *ds1(ds0->GetBlock(blockId));
  vtkDataSet *ds1c(vtkDataSet::SafeDownCast(ds1));
  if(!ds1c)
    throw INTERP_KERNEL::Exception("vtkSedimentDeposit  SplitSingleMultiBloc : nullptr inside single multiblock element !");
  return ds1c;
}

static void ExtractInfo(vtkInformationVector *inputVector, vtkUnstructuredGrid *&usgIn)
{
  vtkInformation *inputInfo(inputVector->GetInformationObject(0));
  vtkDataSet *input = nullptr;
  vtkDataSet *input0(vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  vtkMultiBlockDataSet *input1(vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  if (input0)
  {
    input = input0;
  }
  else
  {
    if (!input1)
    {
      throw INTERP_KERNEL::Exception("Input dataSet must be a DataSet or single elt multi block dataset expected !");
    }
    if (input1->GetNumberOfBlocks() != 1)
    {
      throw INTERP_KERNEL::Exception("Input dataSet is a multiblock dataset with not exactly one block ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    }
    vtkDataObject *input2(input1->GetBlock(0));
    if (!input2)
    {
      throw INTERP_KERNEL::Exception("Input dataSet is a multiblock dataset with exactly one block but this single element is NULL !");
    }
    vtkDataSet *input2c(vtkDataSet::SafeDownCast(input2));
    if (!input2c)
    {
      throw INTERP_KERNEL::Exception("Input dataSet is a multiblock dataset with exactly one block but this single element is not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    }
    input = input2c;
  }
  if (!input)
  {
    throw INTERP_KERNEL::Exception("Input data set is NULL !");
  }
  usgIn = vtkUnstructuredGrid::SafeDownCast(input);
  if (!usgIn)
  {
    throw INTERP_KERNEL::Exception("Input data set is not an unstructured mesh ! This filter works only on unstructured meshes !");
  }
}

std::string vtkSedimentDeposit::vtkInternal::getReprDependingPos(const std::string& origName) const
{
  if( _nb_of_curves == 1 )
    return origName;
  std::ostringstream oss;
  oss << origName << "_" << _curve_id;
  return oss.str();
}

void vtkSedimentDeposit::vtkInternal::fillTable(vtkTable *table) const
{
  if( _curve_id == 0 )
  {
    vtkNew<vtkDoubleArray> timeArr;
    std::string name(getReprDependingPos("Time"));
    timeArr->SetName(name.c_str());
    timeArr->SetNumberOfTuples(_data.size());
    double *pt(timeArr->GetPointer(0));
    {
      std::size_t tmp(0);
      std::for_each(pt, pt + _data.size(), [this, &tmp](double &val) { val = std::get<0>(this->_data[tmp++]); });
    }
    table->AddColumn(timeArr);
  }
  {
    vtkNew<vtkDoubleArray> timeArr;
    std::string name(getReprDependingPos("Total"));
    timeArr->SetName(name.c_str());
    timeArr->SetNumberOfTuples(_data.size());
    double *pt(timeArr->GetPointer(0));
    {
      std::size_t tmp(0);
      std::for_each(pt, pt + _data.size(), [this, &tmp](double &val) { val=std::get<1>(this->_data[tmp])+std::get<2>(this->_data[tmp]); tmp++; });
    }
    table->AddColumn(timeArr);
  }
  {
    vtkNew<vtkDoubleArray> timeArr;
    std::string name(getReprDependingPos("Positif"));
    timeArr->SetName(name.c_str());
    timeArr->SetNumberOfTuples(_data.size());
    double *pt(timeArr->GetPointer(0));
    {
      std::size_t tmp(0);
      std::for_each(pt, pt + _data.size(), [this, &tmp](double &val) { val = std::get<1>(this->_data[tmp++]); });
    }
    table->AddColumn(timeArr);
  }
  {
    vtkNew<vtkDoubleArray> timeArr;
    std::string name(getReprDependingPos("Negatif"));
    timeArr->SetName(name.c_str());
    timeArr->SetNumberOfTuples(_data.size());
    double *pt(timeArr->GetPointer(0));
    {
      std::size_t tmp(0);
      std::for_each(pt, pt + _data.size(), [this, &tmp](double &val) { val = std::get<2>(this->_data[tmp++]); });
    }
    table->AddColumn(timeArr);
  }
}

void vtkSedimentDeposit::vtkInternal::analyzeInputDataSets(vtkUnstructuredGrid *ds1, vtkDataSet *ds2)
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

static std::string Strip(const std::string& ins)
{
  std::string::size_type pos(ins.find_last_not_of(" \t"));
  return ins.substr(0,pos+1);
}

static vtkDataArray *FindFieldWithNameStripped(vtkFieldData *fd, const char *fieldNameToSearch)
{
  std::string keyToSearch(fieldNameToSearch);
  std::vector<std::string> candidates;
  if(!fd)
    throw INTERP_KERNEL::Exception("FindFieldWithNameStripped : nullptr instance !");
  auto nbArrays(fd->GetNumberOfArrays());
  for(auto i = 0 ; i < nbArrays ; ++i)
  {
    std::string arrName(fd->GetArrayName(i));
    if(Strip(arrName) == keyToSearch)
      candidates.push_back(arrName);
  }
  if(candidates.size()!=1)
    {
      std::ostringstream oss; oss << "FindFieldWithNameStripped : not exactly one candidate for \"" << fieldNameToSearch << "\" !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return fd->GetArray(candidates[0].c_str());
}

////////////////////

vtkSedimentDeposit::vtkSedimentDeposit() : NumberOfTimeSteps(0), CurrentTimeIndex(0), IsExecuting(false)
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int vtkSedimentDeposit::RequestUpdateExtent(vtkInformation *info, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
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
  //return vtkDataObjectAlgorithm::RequestUpdateExtent(info,inputVector,outputVector);
}

int vtkSedimentDeposit::RequestInformation(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkSedimentDeposit::RequestInformation ##########################################" << std::endl;
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
    oss << "Exception has been thrown in vtkSedimentDeposit::RequestInformation : " << e.what() << std::endl;
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
  VTKToMEDMemRateOfFlow::WriteMEDFileFromVTKDataSet(mfd, usgIn, {}, 0., 0);
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
  {
    throw INTERP_KERNEL::Exception(message);
  }
}

int vtkSedimentDeposit::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkSedimentDeposit::RequestData        ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid *usgIn(nullptr);
    ExtractInfo(inputVector[0], usgIn);
    vtkInformation *sourceInfo(inputVector[1]->GetInformationObject(0));
    vtkDataObject *source(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
    int nbOfBlocks(GetNumberOfBlocs(source));
    // is this the first request
    if (!this->IsExecuting)
    {
      request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
      this->IsExecuting = true;
      this->Internal2.resize(nbOfBlocks);
      int rk(0);
      std::for_each(this->Internal2.begin(),this->Internal2.end(),[&rk,nbOfBlocks](std::unique_ptr< vtkInternal >& elt) { elt.reset(new vtkInternal(rk++,nbOfBlocks)); });
    }
    //
    this->CurrentTimeIndex++;
    vtkNew<vtkTable> table;
    //
    for(int blockId = 0 ; blockId < nbOfBlocks ; ++blockId)
    {
      vtkDataSet *source1(SplitSingleMultiBloc(source,blockId));
      //
      vtkNew<vtkExplodePolyLine> epl;
      epl->SetInputData(source1);
      epl->Update();
      vtkDataSet *source2(epl->GetOutput());
      //
      this->Internal2[blockId]->analyzeInputDataSets(usgIn, source1);
      ///////////////////////
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> &meshorig(this->Internal2[blockId]->meshOrigin());
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> &untouched_2d_cells(this->Internal2[blockId]->untouched2DCells());
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> &cells_at_boundary_origin(this->Internal2[blockId]->cellsAtBoundary());
      MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> &centers(this->Internal2[blockId]->centers());
      MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> &vol(this->Internal2[blockId]->measure());
      if (this->Internal2[blockId]->computationNeeded())
      {
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> mesh, polygon;
        MEDCoupling::MCAuto<MEDCoupling::MEDFileData> data(MEDCoupling::MEDFileData::New());
        meshorig = ToMedcoupling(data, usgIn);
        meshorig->changeSpaceDimension(2, 0.);
        mesh = meshorig->deepCopy(); //uncouple data and mesh
        {
          MEDCoupling::MCAuto<MEDCoupling::MEDFileData> tmp(MEDCoupling::MEDFileData::New());
          polygon = ToMedcoupling(tmp, source2);
        }
        {
          MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> arr(polygon->getCoords()->keepSelectedComponents({2}));
          if(!arr->isUniform(0, 1e-12))
            std::cerr << "Expected coords array equal to 0 for Z axis... Ignored" << std::endl;
        }
        polygon->changeSpaceDimension(2, 0.);
        {
          bool tmp;
          mcIdType tmpp;
          MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> tmppp(polygon->mergeNodes(1e-12, tmp, tmpp));
        }
        MEDCoupling::MCAuto<MEDCoupling::MEDCoupling1SGTUMesh> polygon_1sgt(MEDCoupling::MEDCoupling1SGTUMesh::New(polygon));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> conn(polygon_1sgt->getNodalConnectivity()->deepCopy());
        conn->rearrange(2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> notNullCellsPolygon;
        {
          MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> conn0(conn->keepSelectedComponents({0})), conn1(conn->keepSelectedComponents({1}));
          MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> delta(MEDCoupling::DataArrayIdType::Substract(conn0, conn1));
          notNullCellsPolygon = delta->findIdsNotEqual(0);
        }
        {
          MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> tmp(polygon_1sgt->buildUnstructured());
          polygon = tmp->buildPartOfMySelf(notNullCellsPolygon->begin(), notNullCellsPolygon->end());
        }
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> polygon_2d(MEDCoupling::MEDCouplingUMesh::New("mesh", 2));
        {
          MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> coo(polygon->getCoords()->deepCopy());
          polygon_2d->setCoords(coo);
        }
        polygon_2d->allocateCells();
        {
          MEDCoupling::MCAuto<MEDCoupling::MEDCoupling1SGTUMesh> tmp(MEDCoupling::MEDCoupling1SGTUMesh::New(polygon));
          conn.takeRef(tmp->getNodalConnectivity());
        }
        conn->rearrange(2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> conn2(conn->fromLinkedListOfPairToList());
        MyAssert(conn2->front() == conn2->back(), "Expected closed wire as input");
        conn2->popBackSilent();
        polygon_2d->insertNextCell(INTERP_KERNEL::NORM_POLYGON, conn2->getNbOfElems(), conn2->begin());
        bool clockWise(false);
        //
        {//false is very important to state if input polygon is clockwise (>0) or not (<0)
          MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> area(polygon_2d->getMeasureField(false));
          clockWise = area->getIJ(0,0) > 0.0;
        }
        //
        double bbox[4];
        mesh->getBoundingBox(bbox);
        double TmpCenter[2] = {(bbox[0] + bbox[1]) / 2., (bbox[2] + bbox[3]) / 2.};
        double Tmpalpha = 1. / std::max(bbox[1] - bbox[0], bbox[3] - bbox[2]);
        double MinusTmpCenter[2] = {-TmpCenter[0], -TmpCenter[1]}, Origin[2] = {0., 0.};
        mesh->translate(MinusTmpCenter);
        mesh->scale(Origin, Tmpalpha);
        polygon->translate(MinusTmpCenter);
        polygon->scale(Origin, Tmpalpha);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> mesh2, line_inter;
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> cellid_in_2d, cellid_in1d;
        {
          MEDCoupling::MEDCouplingUMesh *tmp(nullptr), *tmp2(nullptr);
          MEDCoupling::DataArrayIdType *tmp3(nullptr), *tmp4(nullptr);
          MEDCoupling::MEDCouplingUMesh::Intersect2DMeshWith1DLine(mesh, polygon, 1e-12, tmp, tmp2, tmp3, tmp4);
          mesh2 = tmp;
          line_inter = tmp2;
          cellid_in_2d = tmp3;
          cellid_in1d = tmp4;
        }
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> coo2(mesh2->getCoords()->deepCopy());
        coo2->applyLin(1. / Tmpalpha, 0.);
        coo2->applyLin(1., TmpCenter[0], 0);
        coo2->applyLin(1., TmpCenter[1], 1);
        mesh2->setCoords(coo2);

        std::size_t side = clockWise?1:0;
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> twodcells_to_remove(cellid_in1d->keepSelectedComponents({side}));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> twodcells_to_keep(cellid_in1d->keepSelectedComponents({(side + 1) % 2}));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> ids(twodcells_to_keep->findIdsNotEqual(-1));
        twodcells_to_keep = twodcells_to_keep->selectByTupleId(ids->begin(), ids->end());
        int hotspot(twodcells_to_keep->front());
        twodcells_to_remove = twodcells_to_remove->selectByTupleId(ids->begin(), ids->end());
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> allcells(twodcells_to_remove->buildComplement(mesh2->getNumberOfCells()));
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> mesh2_without_cells_around_polygon(mesh2->buildPartOfMySelf(allcells->begin(), allcells->end()));
        std::vector<MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType>> grps;
        {
          std::vector<MEDCoupling::DataArrayIdType *> tmp(mesh2_without_cells_around_polygon->partitionBySpreadZone());
          std::for_each(tmp.begin(), tmp.end(), [&grps](MEDCoupling::DataArrayIdType *grp) { grps.emplace_back(grp); });
        }
        std::for_each(grps.begin(), grps.end(), [allcells](MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> &grp) { grp = allcells->selectByTupleId(grp->begin(), grp->end()); });
        MyAssert(grps.size() == 2, "Expecting 2 groups, 1 inside, 1 outside !");
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> zeGrp;
        if (grps[0]->presenceOfValue(hotspot) && !grps[1]->presenceOfValue(hotspot))
          zeGrp.takeRef(grps[0]);
        if (grps[1]->presenceOfValue(hotspot) && !grps[0]->presenceOfValue(hotspot))
          zeGrp.takeRef(grps[1]);
        if (zeGrp.isNull())
        {
          throw INTERP_KERNEL::Exception("Internal error : partitioning failed !");
        }
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> mesh3(mesh2->buildPartOfMySelf(zeGrp->begin(), zeGrp->end()));
        double totVol(0.), refVol(0.);
        {
          MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> tmp(mesh3->getMeasureField(true)), tmp2(polygon_2d->getMeasureField(true));
          totVol = tmp->accumulate(0);
          refVol = tmp2->accumulate(0);
        }
        MyAssert(std::abs(totVol - refVol) / refVol < 1e-6, "The test of area conservation failed !");
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> original_cell_ids_2d(cellid_in_2d->selectByTupleId(zeGrp->begin(), zeGrp->end()));
        original_cell_ids_2d->sort(); // les cells dans le referentiel original 2D
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> all_cut_2d_cells(cellid_in1d->deepCopy());
        all_cut_2d_cells->rearrange(1);
        all_cut_2d_cells->sort();
        all_cut_2d_cells = all_cut_2d_cells->buildUnique(); // les cells qui ont subit un split dans le referentiel de sortie 2D
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> all_cut_2d_cells_origin(cellid_in_2d->selectByTupleId(all_cut_2d_cells->begin(), all_cut_2d_cells->end()));
        untouched_2d_cells = original_cell_ids_2d->buildSubstraction(all_cut_2d_cells_origin);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> cells_at_boundary(all_cut_2d_cells->buildIntersection(zeGrp)); // les cellules decoupees dans le referentiel de sortie 2D
        cells_at_boundary_origin = cellid_in_2d->selectByTupleId(cells_at_boundary->begin(), cells_at_boundary->end());
        mesh3 = mesh2->buildPartOfMySelf(cells_at_boundary->begin(), cells_at_boundary->end());
        {
          MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> tmp(mesh3->getMeasureField(true));
          vol.takeRef(tmp->getArray());
        }
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> tmp0(mesh->buildPartOfMySelf(untouched_2d_cells->begin(), untouched_2d_cells->end()));
        mesh3->getBoundingBox(bbox);
        double MinusZeCenter[2] = {-(bbox[0] + bbox[1]) / 2., -(bbox[2] + bbox[3]) / 2.};
        double alpha(1. / std::max(bbox[1] - bbox[0], bbox[3] - bbox[2]));
        mesh3->translate(MinusZeCenter);
        mesh3->scale(Origin, alpha);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> centers_around_zero(mesh3->computeCellCenterOfMass());
        centers = centers_around_zero->deepCopy();
        centers->applyLin(1. / alpha, -MinusZeCenter[0], 0);
        centers->applyLin(1. / alpha, -MinusZeCenter[1], 1);
      }
      constexpr char SEARCHED_FIELD_EVOLUTION[] = "EVOLUTION";
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> f(MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES));
      {
        vtkPointData *pd(usgIn->GetPointData());
        vtkDataArray *evolution_tmp(FindFieldWithNameStripped(pd,SEARCHED_FIELD_EVOLUTION));
        vtkDoubleArray *evolution_tmp2(vtkDoubleArray::SafeDownCast(evolution_tmp));
        if (!evolution_tmp2)
        {
          std::ostringstream oss;
          oss << "Internal error : " << SEARCHED_FIELD_EVOLUTION << " is expected to be of type float32 !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
        MyAssert(evolution_tmp2->GetNumberOfTuples() == meshorig->getNumberOfNodes(), "Mismatch of sizes !");
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> arr(MEDCoupling::DataArrayDouble::New());
        arr->alloc(meshorig->getNumberOfNodes(), 1);
        std::copy(evolution_tmp2->GetPointer(0), evolution_tmp2->GetPointer(meshorig->getNumberOfNodes()), arr->getPointer());
        f->setArray(arr);
        f->setMesh(meshorig);
        f->checkConsistencyLight();
      }
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> f_easy(f->buildSubPart(untouched_2d_cells->begin(), untouched_2d_cells->end()));
      //f_easy->setName("easy");
      //f_easy->writeVTK("easy.vtu");
      MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> weights;
      {
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> tmp(f_easy->getDiscretization()->getMeasureField(f_easy->getMesh(), true));
        weights.takeRef(tmp->getArray());
      }
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> positive_ids(f_easy->getArray()->findIdsGreaterThan(0.));
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> negative_ids(positive_ids->buildComplement(f_easy->getMesh()->getNumberOfNodes()));
      double positive_part(0.), negative_part(0.);
      {
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp(f_easy->getArray()->selectByTupleId(positive_ids->begin(), positive_ids->end()));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp2(weights->selectByTupleId(positive_ids->begin(), positive_ids->end()));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp3(MEDCoupling::DataArrayDouble::Multiply(tmp, tmp2));
        positive_part = tmp3->accumulate((std::size_t)0);
      }
      {
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp(f_easy->getArray()->selectByTupleId(negative_ids->begin(), negative_ids->end()));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp2(weights->selectByTupleId(negative_ids->begin(), negative_ids->end()));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp3(MEDCoupling::DataArrayDouble::Multiply(tmp, tmp2));
        negative_part = tmp3->accumulate((std::size_t)0);
      }
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> f_hard(f->buildSubPart(cells_at_boundary_origin->begin(), cells_at_boundary_origin->end()));
      MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> hard_part(f_hard->getValueOnMulti(centers->begin(), centers->getNumberOfTuples()));
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> positive_ids_hard(hard_part->findIdsGreaterThan(0.));
      MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> negative_ids_hard(positive_ids_hard->buildComplement(cells_at_boundary_origin->getNumberOfTuples()));
      double positive_part_hard(0.), negative_part_hard(0.);
      {
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp(hard_part->selectByTupleId(positive_ids_hard->begin(), positive_ids_hard->end()));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp2(vol->selectByTupleId(positive_ids_hard->begin(), positive_ids_hard->end()));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp3(MEDCoupling::DataArrayDouble::Multiply(tmp, tmp2));
        positive_part_hard = tmp3->accumulate((std::size_t)0);
      }
      {
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp(hard_part->selectByTupleId(negative_ids_hard->begin(), negative_ids_hard->end()));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp2(vol->selectByTupleId(negative_ids_hard->begin(), negative_ids_hard->end()));
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> tmp3(MEDCoupling::DataArrayDouble::Multiply(tmp, tmp2));
        negative_part_hard = tmp3->accumulate((std::size_t)0);
      }
      double timeStep;
      {
        vtkInformation *inInfo(inputVector[0]->GetInformationObject(0));
        vtkDataObject *input(vtkDataObject::GetData(inInfo));
        timeStep = input->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());
      }
      this->Internal2[blockId]->pushData(timeStep, positive_part + positive_part_hard, negative_part + negative_part_hard);
      if (this->CurrentTimeIndex == this->NumberOfTimeSteps)
      {
        this->Internal2[blockId]->fillTable(table);
      }
    }
    if (this->CurrentTimeIndex == this->NumberOfTimeSteps)
      {
        request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
        this->CurrentTimeIndex = 0;
        this->IsExecuting = false;
      }
    vtkInformation *outInfo(outputVector->GetInformationObject(0));
    vtkTable *output(vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
    output->ShallowCopy(table);
  }
  catch (INTERP_KERNEL::Exception &e)
  {
    if (this->IsExecuting)
    {
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      this->CurrentTimeIndex = 0;
      this->IsExecuting = false;
    }
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkSedimentDeposit::RequestData : " << e.what() << std::endl;
    vtkErrorMacro(<< oss.str());
    return 0;
  }
  return 1;
}

void vtkSedimentDeposit::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkSedimentDeposit::SetSourceData(vtkDataObject *input)
{
  this->SetInputData(1, input);
}

void vtkSedimentDeposit::SetSourceConnection(vtkAlgorithmOutput *algOutput)
{
  this->SetInputConnection(1, algOutput);
}

int vtkSedimentDeposit::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
  return 1;
}

bool vtkSedimentDeposit::vtkInternal::computationNeeded() const
{
	if (_recomputationOfMatrixNeeded)
	{
		_meshorig.nullify();
		_untouched_2d_cells.nullify();
		_cells_at_boundary_origin.nullify();
		_centers.nullify();
		_vol.nullify();
	}
	return _recomputationOfMatrixNeeded;
}

