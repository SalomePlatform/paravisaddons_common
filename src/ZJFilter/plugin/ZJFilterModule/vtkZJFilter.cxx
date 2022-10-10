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

#include "vtkZJFilter.h"

#include <vtkAdjacentVertexIterator.h>
#include <vtkAlgorithmOutput.h>
#include <vtkCallbackCommand.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCharArray.h>
#include <vtkMergeBlocks.h>
#include <vtkDataArraySelection.h>
#include <vtkAOSDataArrayTemplate.h>
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
#include <vtkInformationStringKey.h>
#include <vtkInformationVector.h>
#include <vtkLongArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>

#include "InterpKernelException.hxx"
#include "MEDCouplingRefCountObject.hxx"

#include <cstdio>
#include <deque>
#include <map>
#include <sstream>

vtkStandardNewMacro(vtkZJFilter);

///////////////////

vtkInformationDataObjectMetaDataKey* GetMEDReaderMetaDataIfAny()
{
  static const char ZE_KEY[] = "vtkMEDReader::META_DATA";
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

class Grp
{
public:
  Grp(const std::string& name)
    : _name(name)
  {
  }
  void setFamilies(const std::vector<std::string>& fams) { _fams = fams; }
  std::string getName() const { return _name; }
  std::vector<std::string> getFamilies() const { return _fams; }

private:
  std::string _name;
  std::vector<std::string> _fams;
};

class Fam
{
public:
  Fam(const std::string& name)
  {
    constexpr char ZE_SEP[] = "@@][@@";
    std::size_t pos(name.find(ZE_SEP));
    std::string name0(name.substr(0, pos)), name1(name.substr(pos + strlen(ZE_SEP)));
    std::istringstream iss(name1);
    iss >> _id;
    _name = name0;
  }
  std::string getName() const { return _name; }
  int getID() const { return _id; }

private:
  std::string _name;
  int _id;
};

void ExtractInfo(vtkInformationVector* inputVector, vtkUnstructuredGrid*& usgIn)
{
  vtkInformation* inputInfo(inputVector->GetInformationObject(0));
  vtkDataSet* input(0);
  vtkDataSet* input0(vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  vtkMultiBlockDataSet* input1(
    vtkMultiBlockDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  if (input0)
    input = input0;
  else
  {
    if (!input1)
      throw INTERP_KERNEL::Exception(
        "Input dataSet must be a DataSet or single elt multi block dataset expected !");
    if (input1->GetNumberOfBlocks() != 1)
      throw INTERP_KERNEL::Exception(
        "Input dataSet is a multiblock dataset with not exactly one block ! Use MergeBlocks or "
        "ExtractBlocks filter before calling this filter !");
    vtkDataObject* input2(input1->GetBlock(0));
    if (!input2)
      throw INTERP_KERNEL::Exception("Input dataSet is a multiblock dataset with exactly one block "
                                     "but this single element is NULL !");
    vtkDataSet* input2c(vtkDataSet::SafeDownCast(input2));
    if (!input2c)
      throw INTERP_KERNEL::Exception(
        "Input dataSet is a multiblock dataset with exactly one block but this single element is "
        "not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    input = input2c;
  }
  if (!input)
    throw INTERP_KERNEL::Exception("Input data set is NULL !");
  usgIn = vtkUnstructuredGrid::SafeDownCast(input);
  if (!usgIn)
    throw INTERP_KERNEL::Exception("Input data set is not an unstructured mesh ! This filter works "
                                   "only on unstructured meshes !");
}

void LoadFamGrpMapInfo(vtkMutableDirectedGraph* sil, std::string& meshName,
  std::vector<Grp>& groups, std::vector<Fam>& fams)
{
  if (!sil)
    throw INTERP_KERNEL::Exception("LoadFamGrpMapInfo : internal error !");
  int idNames(0);
  vtkAbstractArray* verticesNames(sil->GetVertexData()->GetAbstractArray("Names", idNames));
  vtkStringArray* verticesNames2(vtkStringArray::SafeDownCast(verticesNames));
  vtkIdType id0;
  bool found(false);
  for (int i = 0; i < verticesNames2->GetNumberOfValues(); i++)
  {
    vtkStdString& st(verticesNames2->GetValue(i));
    if (st == "MeshesFamsGrps")
    {
      id0 = i;
      found = true;
    }
  }
  if (!found)
    throw INTERP_KERNEL::Exception(
      "There is an internal error ! The tree on server side has not the expected look !");
  vtkAdjacentVertexIterator* it0(vtkAdjacentVertexIterator::New());
  sil->GetAdjacentVertices(id0, it0);
  int kk(0), ll(0);
  while (it0->HasNext())
  {
    vtkIdType id1(it0->Next());
    std::string mName((const char*)verticesNames2->GetValue(id1));
    meshName = mName;
    vtkAdjacentVertexIterator* it1(vtkAdjacentVertexIterator::New());
    sil->GetAdjacentVertices(id1, it1);
    vtkIdType idZeGrps(it1->Next()); // zeGroups
    vtkAdjacentVertexIterator* itGrps(vtkAdjacentVertexIterator::New());
    sil->GetAdjacentVertices(idZeGrps, itGrps);
    while (itGrps->HasNext())
    {
      vtkIdType idg(itGrps->Next());
      Grp grp((const char*)verticesNames2->GetValue(idg));
      vtkAdjacentVertexIterator* itGrps2(vtkAdjacentVertexIterator::New());
      sil->GetAdjacentVertices(idg, itGrps2);
      std::vector<std::string> famsOnGroup;
      while (itGrps2->HasNext())
      {
        vtkIdType idgf(itGrps2->Next());
        famsOnGroup.push_back(std::string((const char*)verticesNames2->GetValue(idgf)));
      }
      grp.setFamilies(famsOnGroup);
      itGrps2->Delete();
      groups.push_back(grp);
    }
    itGrps->Delete();
    vtkIdType idZeFams(it1->Next()); // zeFams
    it1->Delete();
    vtkAdjacentVertexIterator* itFams(vtkAdjacentVertexIterator::New());
    sil->GetAdjacentVertices(idZeFams, itFams);
    while (itFams->HasNext())
    {
      vtkIdType idf(itFams->Next());
      Fam fam((const char*)verticesNames2->GetValue(idf));
      fams.push_back(fam);
    }
    itFams->Delete();
  }
  it0->Delete();
}

std::vector<std::string> FindConds(const std::vector<Grp>& grps)
{
  constexpr char PAT[] = "COND_";
  constexpr std::size_t SZ_PAT(sizeof(PAT) - 1);
  std::vector<std::string> ret;
  for (std::vector<Grp>::const_iterator it = grps.begin(); it != grps.end(); it++)
  {
    std::string name((*it).getName());
    std::string part(name.substr(0, SZ_PAT));
    if (part == PAT)
      ret.push_back(name.substr(SZ_PAT, std::string::npos));
  }
  return ret;
}

constexpr char EPORT_PAT[] = "EPORT_";

std::vector<std::string> FindEports(
  const std::string& condEntry, const std::vector<Grp>& grps, std::vector<std::string>& eportsZip)
{
  std::vector<std::string> ret;
  std::string commonPart(std::string(EPORT_PAT) + condEntry);
  std::size_t commonPart_sz(commonPart.length());
  for (std::vector<Grp>::const_iterator it = grps.begin(); it != grps.end(); it++)
  {
    std::string name((*it).getName());
    std::string part(name.substr(0, commonPart_sz));
    if (part == commonPart)
    {
      ret.push_back(name);
      eportsZip.push_back(name.substr(commonPart_sz, std::string::npos));
    }
  }
  return ret;
}

std::string BigestCommonPart(const std::string& s1, const std::string& s2)
{
  std::size_t ls1(s1.length()), ls2(s2.length()), lb(0), lt(0);
  std::string b, t;
  if (ls1 >= ls2)
  {
    b = s1;
    t = s2;
    lb = ls1;
    lt = ls2;
  }
  else
  {
    b = s2;
    t = s1;
    lb = ls2;
    lt = ls1;
  }
  for (std::size_t l0 = lt; l0 > 0; l0--)
  {
    for (std::size_t l1 = 0; l1 < lt - l0 + 1; l1++)
    {
      std::string cand(t.substr(l1, l0));
      if (b.find(cand) != std::string::npos)
        return cand;
    }
  }
  return std::string();
}

std::vector<int> DeduceIdsFrom(const std::vector<std::string>& eportsZip)
{
  if (eportsZip.empty())
    return std::vector<int>();
  std::string ref(eportsZip[0]);
  std::size_t sz(eportsZip.size());
  for (std::size_t i = 1; i < sz; i++)
    ref = BigestCommonPart(ref, eportsZip[i]);
  std::vector<int> ret(sz);
  for (std::size_t i = 0; i < sz; i++)
  {
    std::size_t pos(eportsZip[i].find(ref));
    if (pos == std::string::npos)
      throw INTERP_KERNEL::Exception("DeduceIdsFrom : internal error !");
    std::string res(
      eportsZip[i].substr(0, pos) + eportsZip[i].substr(pos + ref.length(), std::string::npos));
    std::istringstream iss(res);
    int val(0);
    iss >> val;
    ret[i] = val;
  }
  return ret;
}

std::set<int> FamiliesIdsFromGrp(
  const std::vector<Grp>& grps, const std::vector<Fam>& fams, const std::string& grp)
{
  bool found(false);
  std::vector<std::string> locFams;
  for (std::vector<Grp>::const_iterator it = grps.begin(); it != grps.end(); it++)
  {
    if ((*it).getName() == grp)
    {
      locFams = (*it).getFamilies();
      found = true;
      break;
    }
  }
  if (!found)
    throw INTERP_KERNEL::Exception("FamiliesIdsFromGrp : internal error !");
  std::set<int> ret;
  for (std::vector<Fam>::const_iterator it = fams.begin(); it != fams.end(); it++)
  {
    if (std::find(locFams.begin(), locFams.end(), (*it).getName()) != locFams.end())
      ret.insert((*it).getID());
  }
  return ret;
}

vtkDataSet* FilterFamilies(vtkZJFilter* zeBoss, vtkDataSet* input, const std::set<int>& idsToKeep)
{
  bool catchAll, catchSmth;
  vtkNew<vtkThreshold> thres;
  thres->AddObserver(vtkCommand::ProgressEvent, zeBoss->InternalProgressObserver);
  constexpr int VTK_DATA_ARRAY_DELETE = vtkAOSDataArrayTemplate<double>::VTK_DATA_ARRAY_DELETE;
  constexpr char ZE_SELECTION_ARR_NAME[] = "@@ZeSelection@@";
  constexpr char arrNameOfFamilyField[] = "FamilyIdCell";
  constexpr char associationForThreshold[] = "vtkDataObject::FIELD_ASSOCIATION_CELLS";
  vtkDataSet* output(input->NewInstance());
  output->ShallowCopy(input);
  thres->SetInputData(output);
  vtkDataSetAttributes *dscIn(input->GetCellData()), *dscIn2(input->GetPointData());
  vtkDataSetAttributes *dscOut(output->GetCellData()), *dscOut2(output->GetPointData());
  //
  constexpr double vMin(1.), vMax(2.);
  thres->SetUpperThreshold(vMax);
  thres->SetLowerThreshold(vMin);
  // OK for the output
  //
  vtkDataArray* da(input->GetCellData()->GetScalars(arrNameOfFamilyField));
  if (!da)
    return 0;
  std::string daName(da->GetName());
  vtkLongArray* dai(vtkLongArray::SafeDownCast(da));
  if (daName != arrNameOfFamilyField || !dai)
    return 0;
  //
  int nbOfTuples(dai->GetNumberOfTuples());
  vtkCharArray* zeSelection(vtkCharArray::New());
  zeSelection->SetName(ZE_SELECTION_ARR_NAME);
  zeSelection->SetNumberOfComponents(1);
  char* pt(new char[nbOfTuples]);
  zeSelection->SetArray(pt, nbOfTuples, 0, VTK_DATA_ARRAY_DELETE);
  const long* inPtr(dai->GetPointer(0));
  std::fill(pt, pt + nbOfTuples, 0);
  catchAll = true;
  catchSmth = false;
  std::vector<bool> pt2(nbOfTuples, false);
  for (std::set<int>::const_iterator it = idsToKeep.begin(); it != idsToKeep.end(); it++)
  {
    bool catchFid(false);
    for (int i = 0; i < nbOfTuples; i++)
      if (inPtr[i] == *it)
      {
        pt2[i] = true;
        catchFid = true;
      }
    if (!catchFid)
      catchAll = false;
    else
      catchSmth = true;
  }
  for (int ii = 0; ii < nbOfTuples; ii++)
    if (pt2[ii])
      pt[ii] = 2;
  int idx(output->GetCellData()->AddArray(zeSelection));
  output->GetCellData()->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  output->GetCellData()->CopyScalarsOff();
  zeSelection->Delete();
  //
  thres->SetInputArrayToProcess(idx, 0, 0, associationForThreshold, ZE_SELECTION_ARR_NAME);
  thres->Update();
  vtkUnstructuredGrid* zeComputedOutput(thres->GetOutput());
  zeComputedOutput->GetCellData()->RemoveArray(idx);
  output->Delete();
  zeComputedOutput->Register(0);
  thres->RemoveObserver(zeBoss->InternalProgressObserver);
  return zeComputedOutput;
}

////////////////////

vtkZJFilter::vtkZJFilter()
  : InternalProgressObserver(0)
{
  this->InternalProgressObserver = vtkCallbackCommand::New();
  this->InternalProgressObserver->SetCallback(&vtkZJFilter::InternalProgressCallbackFunction);
  this->InternalProgressObserver->SetClientData(this);
}

vtkZJFilter::~vtkZJFilter()
{
  this->InternalProgressObserver->Delete();
}

void vtkZJFilter::InternalProgressCallbackFunction(
  vtkObject* arg, unsigned long, void* clientdata, void*)
{
  reinterpret_cast<vtkZJFilter*>(clientdata)
    ->InternalProgressCallback(static_cast<vtkAlgorithm*>(arg));
}

void vtkZJFilter::InternalProgressCallback(vtkAlgorithm* algorithm)
{
  /*this->UpdateProgress(algorithm->GetProgress()); // To intercept progression of Threshold filters
  if (this->AbortExecute)
    {
      algorithm->SetAbortExecute(1);
      }*/
}

int vtkZJFilter::RequestInformation(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // std::cerr << "########################################## vtkZJFilter::RequestInformation
  // ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid* usgIn = nullptr;
    ExtractInfo(inputVector[0], usgIn);
  }
  catch (INTERP_KERNEL::Exception& e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkZJFilter::RequestInformation : " << e.what()
        << std::endl;
    if (this->HasObserver("ErrorEvent"))
    {
      this->InvokeEvent("ErrorEvent", const_cast<char*>(oss.str().c_str()));
    }
    else
    {
      vtkOutputWindowDisplayErrorText(const_cast<char*>(oss.str().c_str()));
    }
    vtkObject::BreakOnError();
    return 0;
  }
  return 1;
}

int vtkZJFilter::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // std::cerr << "########################################## vtkZJFilter::RequestData
  // ##########################################" << std::endl;
  try
  {
    vtkInformation *inputInfo(inputVector[0]->GetInformationObject(0));
    vtkInformation *outInfo(outputVector->GetInformationObject(0));
    vtkUnstructuredGrid* usgIn(nullptr);
    ExtractInfo(inputVector[0], usgIn);
    std::string meshName;
    std::vector<Grp> groups;
    std::vector<Fam> fams;
    if (IsInformationOK(inputInfo))
    {
      vtkMutableDirectedGraph* famGrpGraph(
        vtkMutableDirectedGraph::SafeDownCast(inputInfo->Get(GetMEDReaderMetaDataIfAny())));
      LoadFamGrpMapInfo(famGrpGraph, meshName, groups, fams);
    }
    std::vector<std::string> conds(FindConds(groups));
    vtkUnstructuredGrid* output(
      vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
    vtkNew<vtkMultiBlockDataSet> mb2;
    std::size_t iblock(0);
    for (std::vector<std::string>::const_iterator it = conds.begin(); it != conds.end();
         it++, iblock++)
    {
      std::vector<std::string> eports2;
      std::vector<std::string> eports(FindEports(*it, groups, eports2));
      std::vector<int> ids(DeduceIdsFrom(eports2));
      vtkNew<vtkMultiBlockDataSet> mb;
      std::size_t sz(eports.size());
      for (std::size_t i = 0; i < sz; i++)
      {
        std::set<int> zeIds(FamiliesIdsFromGrp(groups, fams, eports[i]));
        //
        vtkSmartPointer<vtkDataSet> ds(FilterFamilies(this, usgIn, zeIds));
        {
          vtkNew<vtkLongArray> arr;
          arr->SetName((*it).c_str());
          arr->SetNumberOfComponents(1);
          int nbTuples(ds->GetNumberOfCells());
          arr->SetNumberOfTuples(nbTuples);
          long* pt(arr->GetPointer(0));
          std::fill(pt, pt + nbTuples, ids[i]);
          ds->GetCellData()->AddArray(arr);
        }
        this->UpdateProgress(double(i) / double(sz));
        mb->SetBlock(i, ds);
      }
      vtkNew<vtkMergeBlocks> cd;
      cd->SetInputData(mb);
      cd->SetMergePoints(0);
      cd->Update();
      mb2->SetBlock(iblock, cd->GetOutput());
    }
    {
      vtkNew<vtkMergeBlocks> cd2;
      cd2->SetInputData(mb2);
      cd2->SetMergePoints(0);
      cd2->Update();
      output->ShallowCopy(cd2->GetOutput());
    }
  }
  catch (INTERP_KERNEL::Exception& e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkZJFilter::RequestInformation : " << e.what()
        << std::endl;
    if (this->HasObserver("ErrorEvent"))
      this->InvokeEvent("ErrorEvent", const_cast<char*>(oss.str().c_str()));
    else
      vtkOutputWindowDisplayErrorText(const_cast<char*>(oss.str().c_str()));
    vtkObject::BreakOnError();
    return 0;
  }
  return 1;
}

void vtkZJFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
