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
// Author : Anthony Geay

#include "vtkExtractThreeD.h"

#include <NormalizedGeometricTypes>
#include <CellModel.hxx>
#include <InterpKernelException.hxx>
#include <MEDFileFieldOverView.hxx>

#include <vtkAdjacentVertexIterator.h>
#include <vtkAlgorithmOutput.h>
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArraySelection.h>
#include <vtkAOSDataArrayTemplate.h>
#include <vtkDataObjectTreeIterator.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkInEdgeIterator.h>
#include <vtkInformation.h>
#include <vtkInformationDataObjectKey.h>
#include <vtkInformationStringKey.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkThreshold.h>
#include <vtkTimeStamp.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVariantArray.h>

#include <deque>
#include <map>
#include <sstream>

vtkStandardNewMacro(vtkExtractThreeD);

///////////////////

class vtkExtractThreeD::vtkExtractThreeDInternal
{
public:
  vtkExtractThreeDInternal() {}
  std::vector<int> getIdsToKeep() const { return _types; }
  void setIdsToKeep(const std::vector<int> &types) { _types = types; }

private:
  std::vector<int> _types;
};

vtkExtractThreeD::vtkExtractThreeD() : Internal(new vtkExtractThreeDInternal)
{
}

vtkExtractThreeD::~vtkExtractThreeD()
{
  delete this->Internal;
}

int vtkExtractThreeD::RequestInformation(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  try
  {
    //std::cerr << "########################################## vtkExtractThreeD::RequestInformation ##########################################" << std::endl;
    vtkInformation *outInfo(outputVector->GetInformationObject(0));
    vtkInformation *inputInfo(inputVector[0]->GetInformationObject(0));
    vtkDataSet *input(0);
    {
      vtkDataObject *inp(inputInfo->Get(vtkDataObject::DATA_OBJECT()));
      if (vtkDataSet::SafeDownCast(inp))
        input = vtkDataSet::SafeDownCast(inp);
      else
      {
        vtkMultiBlockDataSet *inputTmp(vtkMultiBlockDataSet::SafeDownCast(inp));
        if (inputTmp)
        {
          if (inputTmp->GetNumberOfBlocks() != 1)
          {
            vtkDebugMacro("vtkExtractThreeD::RequestInformation : input vtkMultiBlockDataSet must contain exactly 1 block !");
            return 0;
          }
          vtkDataSet *blk0(vtkDataSet::SafeDownCast(inputTmp->GetBlock(0)));
          if (!blk0)
          {
            vtkDebugMacro("vtkExtractThreeD::RequestInformation : the single block in input vtkMultiBlockDataSet must be a vtkDataSet instance !");
            return 0;
          }
          input = blk0;
        }
        else
        {
          vtkDebugMacro("vtkExtractThreeD::RequestInformation : supported input are vtkDataSet or vtkMultiBlockDataSet !");
          return 0;
        }
      }
    }
    {
      vtkIdType nbOfCells(input->GetNumberOfCells());
      std::map<int, INTERP_KERNEL::NormalizedCellType> m;
      std::set<int> typesToKeep;
      for (vtkIdType cellId = 0; cellId < nbOfCells; cellId++)
      {
        int vtkCt(input->GetCellType(cellId));
        const std::map<int, INTERP_KERNEL::NormalizedCellType>::const_iterator it(m.find(vtkCt));
        if (it == m.end())
        {
          const unsigned char *pos(std::find(MEDCoupling::MEDMeshMultiLev::PARAMEDMEM_2_VTKTYPE, MEDCoupling::MEDMeshMultiLev::PARAMEDMEM_2_VTKTYPE + MEDCoupling::MEDMeshMultiLev::PARAMEDMEM_2_VTKTYPE_LGTH, vtkCt));
          if (pos == MEDCoupling::MEDMeshMultiLev::PARAMEDMEM_2_VTKTYPE + MEDCoupling::MEDMeshMultiLev::PARAMEDMEM_2_VTKTYPE_LGTH)
          {
            vtkDebugMacro("vtkExtractThreeD::RequestInformation : cell #" << cellId << " has unrecognized type !");
            return 0;
          }
          INTERP_KERNEL::NormalizedCellType mcCtype((INTERP_KERNEL::NormalizedCellType)std::distance(MEDCoupling::MEDMeshMultiLev::PARAMEDMEM_2_VTKTYPE, pos));
          const INTERP_KERNEL::CellModel &cm(INTERP_KERNEL::CellModel::GetCellModel(mcCtype));
          if (cm.getDimension() == 3)
            typesToKeep.insert(vtkCt);
        }
      }
      std::vector<int> typesToKeep2(typesToKeep.begin(), typesToKeep.end());
      this->Internal->setIdsToKeep(typesToKeep2);
    }
  }
  catch (INTERP_KERNEL::Exception &e)
  {
    std::cerr << "Exception has been thrown in vtkExtractThreeD::RequestInformation : " << e.what() << std::endl;
    return 0;
  }
  return 1;
}

vtkDataSet *FilterFamilies(vtkDataSet *input, const std::vector<int> &idsToKeep)
{
  const int VTK_DATA_ARRAY_DELETE = vtkAOSDataArrayTemplate<double>::VTK_DATA_ARRAY_DELETE;
  const char ZE_SELECTION_ARR_NAME[] = "@@ZeSelection@@";
  vtkDataSet *output(input->NewInstance());
  output->ShallowCopy(input);
  vtkSmartPointer<vtkThreshold> thres(vtkSmartPointer<vtkThreshold>::New());
  thres->SetInputData(output);
  vtkDataSetAttributes *dscIn(input->GetCellData()), *dscIn2(input->GetPointData());
  vtkDataSetAttributes *dscOut(output->GetCellData()), *dscOut2(output->GetPointData());
  //
  double vMin(1.), vMax(2.);
  thres->SetUpperThreshold(vMax);
  thres->SetLowerThreshold(vMin);
  // OK for the output
  vtkIdType nbOfCells(input->GetNumberOfCells());
  vtkCharArray *zeSelection(vtkCharArray::New());
  zeSelection->SetName(ZE_SELECTION_ARR_NAME);
  zeSelection->SetNumberOfComponents(1);
  char *pt(new char[nbOfCells]);
  zeSelection->SetArray(pt, nbOfCells, 0, VTK_DATA_ARRAY_DELETE);
  std::fill(pt, pt + nbOfCells, 0);
  std::vector<bool> pt2(nbOfCells, false);
  for (std::vector<int>::const_iterator it = idsToKeep.begin(); it != idsToKeep.end(); it++)
  {
    for (vtkIdType ii = 0; ii < nbOfCells; ii++)
    {
      if (input->GetCellType(ii) == *it)
        pt2[ii] = true;
    }
  }
  for (int ii = 0; ii < nbOfCells; ii++)
    if (pt2[ii])
      pt[ii] = 2;
  int idx(output->GetCellData()->AddArray(zeSelection));
  output->GetCellData()->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  output->GetCellData()->CopyScalarsOff();
  zeSelection->Delete();
  //
  thres->SetInputArrayToProcess(idx, 0, 0, "vtkDataObject::FIELD_ASSOCIATION_CELLS", ZE_SELECTION_ARR_NAME);
  thres->Update();
  vtkUnstructuredGrid *zeComputedOutput(thres->GetOutput());
  zeComputedOutput->GetCellData()->RemoveArray(idx);
  output->Delete();
  zeComputedOutput->Register(0);
  return zeComputedOutput;
}

int vtkExtractThreeD::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  try
  {
    //std::cerr << "########################################## vtkExtractThreeD::RequestData        ##########################################" << std::endl;
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);
    vtkDataSet *input(vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
    vtkInformation *info(input->GetInformation());
    vtkInformation *outInfo(outputVector->GetInformationObject(0));
    vtkDataSet *output(vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
    std::vector<int> idsToKeep(this->Internal->getIdsToKeep());
    vtkDataSet *tryOnCell(FilterFamilies(input, idsToKeep));
    // first shrink the input
    output->ShallowCopy(tryOnCell);
    tryOnCell->Delete();
  }
  catch (INTERP_KERNEL::Exception &e)
  {
    std::cerr << "Exception has been thrown in vtkExtractThreeD::RequestData : " << e.what() << std::endl;
    return 0;
  }
  return 1;
}

void vtkExtractThreeD::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
