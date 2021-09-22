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

#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelGaussCoords.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingUMesh.hxx"

#include <deque>
#include <map>
#include <set>
#include <sstream>

using MEDCoupling::DataArray;
using MEDCoupling::DataArrayDouble;
using MEDCoupling::DataArrayInt;
using MEDCoupling::DataArrayInt64;
using MEDCoupling::DynamicCastSafe;
using MEDCoupling::MCAuto;
using MEDCoupling::MEDCouplingFieldDouble;
using MEDCoupling::MEDCouplingMesh;
using MEDCoupling::MEDCouplingUMesh;
using MEDCoupling::ON_GAUSS_PT;

vtkStandardNewMacro(vtkTorseurCIH);
///////////////////

std::map<int, int> ComputeMapOfType()
{
  std::map<int, int> ret;
  int nbOfTypesInMC(sizeof(MEDCOUPLING2VTKTYPETRADUCER) / sizeof(int));
  for (int i = 0; i < nbOfTypesInMC; i++)
  {
    int vtkId(MEDCOUPLING2VTKTYPETRADUCER[i]);
    if (vtkId != -1)
      ret[vtkId] = i;
  }
  return ret;
}

std::map<int, int> ComputeRevMapOfType()
{
  std::map<int, int> ret;
  int nbOfTypesInMC(sizeof(MEDCOUPLING2VTKTYPETRADUCER) / sizeof(int));
  for (int i = 0; i < nbOfTypesInMC; i++)
  {
    int vtkId(MEDCOUPLING2VTKTYPETRADUCER[i]);
    if (vtkId != -1)
      ret[i] = vtkId;
  }
  return ret;
}

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

DataArrayInt* ConvertVTKArrayToMCArrayInt(vtkDataArray* data)
{
  if (!data)
    throw INTERP_KERNEL::Exception("ConvertVTKArrayToMCArrayInt : internal error !");
  int nbTuples(data->GetNumberOfTuples()), nbComp(data->GetNumberOfComponents());
  std::size_t nbElts(nbTuples * nbComp);
  MCAuto<DataArrayInt> ret(DataArrayInt::New());
  ret->alloc(nbTuples, nbComp);
  for (int i = 0; i < nbComp; i++)
  {
    const char* comp(data->GetComponentName(i));
    if (comp)
      ret->setInfoOnComponent(i, comp);
  }
  int* ptOut(ret->getPointer());
  vtkIntArray* d0(vtkIntArray::SafeDownCast(data));
  if (d0)
  {
    const int* pt(d0->GetPointer(0));
    std::copy(pt, pt + nbElts, ptOut);
    return ret.retn();
  }
  vtkLongArray* d1(vtkLongArray::SafeDownCast(data));
  if (d1)
  {
    const long* pt(d1->GetPointer(0));
    std::copy(pt, pt + nbElts, ptOut);
    return ret.retn();
  }
  vtkIdTypeArray* d2(vtkIdTypeArray::SafeDownCast(data));
  if (d2)
  {
    const vtkIdType* pt(d2->GetPointer(0));
    std::copy(pt, pt + nbElts, ptOut);
    return ret.retn();
  }
  std::ostringstream oss;
  oss << "ConvertVTKArrayToMCArrayInt : unrecognized array \"" << typeid(*data).name()
      << "\" type !";
  throw INTERP_KERNEL::Exception(oss.str());
}

DataArrayDouble* ConvertVTKArrayToMCArrayDouble(vtkDataArray* data)
{
  if (!data)
    throw INTERP_KERNEL::Exception("ConvertVTKArrayToMCArrayDouble : internal error !");
  int nbTuples(data->GetNumberOfTuples()), nbComp(data->GetNumberOfComponents());
  std::size_t nbElts(nbTuples * nbComp);
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  ret->alloc(nbTuples, nbComp);
  for (int i = 0; i < nbComp; i++)
  {
    const char* comp(data->GetComponentName(i));
    if (comp)
      ret->setInfoOnComponent(i, comp);
  }
  double* ptOut(ret->getPointer());
  vtkFloatArray* d0(vtkFloatArray::SafeDownCast(data));
  if (d0)
  {
    const float* pt(d0->GetPointer(0));
    for (std::size_t i = 0; i < nbElts; i++)
      ptOut[i] = pt[i];
    return ret.retn();
  }
  vtkDoubleArray* d1(vtkDoubleArray::SafeDownCast(data));
  if (d1)
  {
    const double* pt(d1->GetPointer(0));
    std::copy(pt, pt + nbElts, ptOut);
    return ret.retn();
  }
  std::ostringstream oss;
  oss << "ConvertVTKArrayToMCArrayDouble : unrecognized array \"" << typeid(*data).name()
      << "\" type !";
  throw INTERP_KERNEL::Exception(oss.str());
}

DataArray* ConvertVTKArrayToMCArray(vtkDataArray* data)
{
  if (!data)
    throw INTERP_KERNEL::Exception("ConvertVTKArrayToMCArray : internal error !");
  vtkFloatArray* d0(vtkFloatArray::SafeDownCast(data));
  vtkDoubleArray* d1(vtkDoubleArray::SafeDownCast(data));
  if (d0 || d1)
    return ConvertVTKArrayToMCArrayDouble(data);
  vtkIntArray* d2(vtkIntArray::SafeDownCast(data));
  vtkLongArray* d3(vtkLongArray::SafeDownCast(data));
  if (d2 || d3)
    return ConvertVTKArrayToMCArrayInt(data);
  std::ostringstream oss;
  oss << "ConvertVTKArrayToMCArray : unrecognized array \"" << typeid(*data).name() << "\" type !";
  throw INTERP_KERNEL::Exception(oss.str());
}

DataArrayDouble* BuildCoordsFrom(vtkPointSet* ds)
{
  if (!ds)
    throw INTERP_KERNEL::Exception("BuildCoordsFrom : internal error !");
  vtkPoints* pts(ds->GetPoints());
  if (!pts)
    throw INTERP_KERNEL::Exception("BuildCoordsFrom : internal error 2 !");
  vtkDataArray* data(pts->GetData());
  if (!data)
    throw INTERP_KERNEL::Exception("BuildCoordsFrom : internal error 3 !");
  MCAuto<DataArrayDouble> coords(ConvertVTKArrayToMCArrayDouble(data));
  return coords.retn();
}

void ConvertFromUnstructuredGrid(vtkUnstructuredGrid* ds,
  std::vector<MCAuto<MEDCouplingUMesh> >& ms, std::vector<MCAuto<DataArrayIdType> >& ids)
{
  MCAuto<DataArrayDouble> coords(BuildCoordsFrom(ds));
  vtkIdType nbCells(ds->GetNumberOfCells());
  vtkUnsignedCharArray* ct(ds->GetCellTypesArray());
  if (!ct)
    throw INTERP_KERNEL::Exception("ConvertFromUnstructuredGrid : internal error");
  const unsigned char* ctPtr(ct->GetPointer(0));
  std::map<int, int> m(ComputeMapOfType());
  MCAuto<DataArrayIdType> lev(DataArrayIdType::New());
  lev->alloc(nbCells, 1);
  mcIdType* levPtr(lev->getPointer());
  for (vtkIdType i = 0; i < nbCells; i++)
  {
    std::map<int, int>::iterator it(m.find(ctPtr[i]));
    if (it != m.end())
    {
      const INTERP_KERNEL::CellModel& cm(
        INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)(*it).second));
      levPtr[i] = cm.getDimension();
    }
    else
    {
      std::ostringstream oss;
      oss << "ConvertFromUnstructuredGrid : at pos #" << i
          << " unrecognized VTK cell with type =" << ctPtr[i];
      throw INTERP_KERNEL::Exception(oss.str());
    }
  }
  MCAuto<DataArrayIdType> levs(lev->getDifferentValues());
  vtkIdTypeArray *faces(ds->GetFaces()), *faceLoc(ds->GetFaceLocations());
  for (const mcIdType* curLev = levs->begin(); curLev != levs->end(); curLev++)
  {
    MCAuto<MEDCouplingUMesh> m0(MEDCouplingUMesh::New("", *curLev));
    m0->setCoords(coords);
    m0->allocateCells();
    MCAuto<DataArrayIdType> cellIdsCurLev(lev->findIdsEqual(*curLev));
    for (const mcIdType* cellId = cellIdsCurLev->begin(); cellId != cellIdsCurLev->end(); cellId++)
    {
      std::map<int, int>::iterator it(m.find(ctPtr[*cellId]));
      vtkIdType sz;
      vtkIdType const* pts;
      ds->GetCellPoints(*cellId, sz, pts);
      INTERP_KERNEL::NormalizedCellType ct((INTERP_KERNEL::NormalizedCellType)(*it).second);
      if (ct != INTERP_KERNEL::NORM_POLYHED)
      {
        std::vector<mcIdType> conn2(sz);
        for (int kk = 0; kk < sz; kk++)
          conn2[kk] = pts[kk];
        m0->insertNextCell(ct, sz, &conn2[0]);
      }
      else
      {
        if (!faces || !faceLoc)
          throw INTERP_KERNEL::Exception(
            "ConvertFromUnstructuredGrid : faces are expected when there are polyhedra !");
        const vtkIdType *facPtr(faces->GetPointer(0)), *facLocPtr(faceLoc->GetPointer(0));
        std::vector<mcIdType> conn;
        int off0(facLocPtr[*cellId]);
        int nbOfFaces(facPtr[off0++]);
        for (int k = 0; k < nbOfFaces; k++)
        {
          int nbOfNodesInFace(facPtr[off0++]);
          std::copy(facPtr + off0, facPtr + off0 + nbOfNodesInFace, std::back_inserter(conn));
          off0 += nbOfNodesInFace;
          if (k < nbOfFaces - 1)
            conn.push_back(-1);
        }
        m0->insertNextCell(ct, conn.size(), conn.data());
      }
    }
    ms.push_back(m0);
    ids.push_back(cellIdsCurLev);
  }
}

vtkSmartPointer<vtkUnstructuredGrid> ConvertUMeshFromMCToVTK(const MEDCouplingUMesh* mVor)
{
  std::map<int, int> zeMapRev(ComputeRevMapOfType());
  int nbCells(mVor->getNumberOfCells());
  vtkSmartPointer<vtkUnstructuredGrid> ret(vtkSmartPointer<vtkUnstructuredGrid>::New());
  ret->Initialize();
  ret->Allocate();
  vtkSmartPointer<vtkPoints> points(vtkSmartPointer<vtkPoints>::New());
  {
    const DataArrayDouble* vorCoords(mVor->getCoords());
    vtkSmartPointer<vtkDoubleArray> da(vtkSmartPointer<vtkDoubleArray>::New());
    da->SetNumberOfComponents(vorCoords->getNumberOfComponents());
    da->SetNumberOfTuples(vorCoords->getNumberOfTuples());
    std::copy(vorCoords->begin(), vorCoords->end(), da->GetPointer(0));
    points->SetData(da);
  }
  mVor->checkConsistencyLight();
  switch (mVor->getMeshDimension())
  {
    case 3:
    {
      vtkIdType *cPtr(nullptr), *dPtr(nullptr);
      unsigned char* aPtr(nullptr);
      vtkSmartPointer<vtkUnsignedCharArray> cellTypes(vtkSmartPointer<vtkUnsignedCharArray>::New());
      {
        cellTypes->SetNumberOfComponents(1);
        cellTypes->SetNumberOfTuples(nbCells);
        aPtr = cellTypes->GetPointer(0);
      }
      vtkSmartPointer<vtkIdTypeArray> cellLocations(vtkSmartPointer<vtkIdTypeArray>::New());
      {
        cellLocations->SetNumberOfComponents(1);
        cellLocations->SetNumberOfTuples(nbCells);
        cPtr = cellLocations->GetPointer(0);
      }
      vtkSmartPointer<vtkIdTypeArray> cells(vtkSmartPointer<vtkIdTypeArray>::New());
      {
        MCAuto<DataArrayIdType> tmp2(mVor->computeEffectiveNbOfNodesPerCell());
        cells->SetNumberOfComponents(1);
        cells->SetNumberOfTuples(tmp2->accumulate((std::size_t)0) + nbCells);
        dPtr = cells->GetPointer(0);
      }
      const mcIdType *connPtr(mVor->getNodalConnectivity()->begin()),
        *connIPtr(mVor->getNodalConnectivityIndex()->begin());
      int k(0), kk(0);
      std::vector<vtkIdType> ee, ff;
      for (int i = 0; i < nbCells; i++, connIPtr++)
      {
        INTERP_KERNEL::NormalizedCellType ct(
          static_cast<INTERP_KERNEL::NormalizedCellType>(connPtr[connIPtr[0]]));
        *aPtr++ = zeMapRev[connPtr[connIPtr[0]]];
        if (ct != INTERP_KERNEL::NORM_POLYHED)
        {
          int sz(connIPtr[1] - connIPtr[0] - 1);
          *dPtr++ = sz;
          dPtr = std::copy(connPtr + connIPtr[0] + 1, connPtr + connIPtr[1], dPtr);
          *cPtr++ = k;
          k += sz + 1;
          ee.push_back(kk);
        }
        else
        {
          std::set<mcIdType> s(connPtr + connIPtr[0] + 1, connPtr + connIPtr[1]);
          s.erase(-1);
          int nbFace(std::count(connPtr + connIPtr[0] + 1, connPtr + connIPtr[1], -1) + 1);
          ff.push_back(nbFace);
          const mcIdType* work(connPtr + connIPtr[0] + 1);
          for (int j = 0; j < nbFace; j++)
          {
            const mcIdType* work2 = std::find(work, connPtr + connIPtr[1], -1);
            ff.push_back(std::distance(work, work2));
            ff.insert(ff.end(), work, work2);
            work = work2 + 1;
          }
          *dPtr++ = (int)s.size();
          dPtr = std::copy(s.begin(), s.end(), dPtr);
          *cPtr++ = k;
          k += (int)s.size() + 1;
          ee.push_back(kk);
          kk += connIPtr[1] - connIPtr[0] + 1;
        }
      }
      //
      vtkSmartPointer<vtkIdTypeArray> faceLocations(vtkSmartPointer<vtkIdTypeArray>::New());
      {
        faceLocations->SetNumberOfComponents(1);
        faceLocations->SetNumberOfTuples(ee.size());
        std::copy(ee.begin(), ee.end(), faceLocations->GetPointer(0));
      }
      vtkSmartPointer<vtkIdTypeArray> faces(vtkSmartPointer<vtkIdTypeArray>::New());
      {
        faces->SetNumberOfComponents(1);
        faces->SetNumberOfTuples(ff.size());
        std::copy(ff.begin(), ff.end(), faces->GetPointer(0));
      }
      vtkSmartPointer<vtkCellArray> cells2(vtkSmartPointer<vtkCellArray>::New());
      cells2->SetCells(nbCells, cells);
      ret->SetCells(cellTypes, cellLocations, cells2, faceLocations, faces);
      break;
    }
    case 2:
    {
      vtkSmartPointer<vtkUnsignedCharArray> cellTypes(vtkSmartPointer<vtkUnsignedCharArray>::New());
      {
        cellTypes->SetNumberOfComponents(1);
        cellTypes->SetNumberOfTuples(nbCells);
        unsigned char* ptr(cellTypes->GetPointer(0));
        std::fill(ptr, ptr + nbCells, zeMapRev[(int)INTERP_KERNEL::NORM_POLYGON]);
      }
      vtkIdType *cPtr(nullptr), *dPtr(nullptr);
      vtkSmartPointer<vtkIdTypeArray> cellLocations(vtkSmartPointer<vtkIdTypeArray>::New());
      {
        cellLocations->SetNumberOfComponents(1);
        cellLocations->SetNumberOfTuples(nbCells);
        cPtr = cellLocations->GetPointer(0);
      }
      vtkSmartPointer<vtkIdTypeArray> cells(vtkSmartPointer<vtkIdTypeArray>::New());
      {
        cells->SetNumberOfComponents(1);
        cells->SetNumberOfTuples(mVor->getNodalConnectivity()->getNumberOfTuples());
        dPtr = cells->GetPointer(0);
      }
      const mcIdType *connPtr(mVor->getNodalConnectivity()->begin()),
        *connIPtr(mVor->getNodalConnectivityIndex()->begin());
      mcIdType k(0);
      for (mcIdType i = 0; i < nbCells; i++, connIPtr++)
      {
        *dPtr++ = connIPtr[1] - connIPtr[0] - 1;
        dPtr = std::copy(connPtr + connIPtr[0] + 1, connPtr + connIPtr[1], dPtr);
        *cPtr++ = k;
        k += connIPtr[1] - connIPtr[0];
      }
      vtkSmartPointer<vtkCellArray> cells2(vtkSmartPointer<vtkCellArray>::New());
      cells2->SetCells(nbCells, cells);
      ret->SetCells(cellTypes, cellLocations, cells2);
      break;
    }
    case 1:
    {
      vtkSmartPointer<vtkUnsignedCharArray> cellTypes(vtkSmartPointer<vtkUnsignedCharArray>::New());
      {
        cellTypes->SetNumberOfComponents(1);
        cellTypes->SetNumberOfTuples(nbCells);
        unsigned char* ptr(cellTypes->GetPointer(0));
        std::fill(ptr, ptr + nbCells, zeMapRev[(int)INTERP_KERNEL::NORM_SEG2]);
      }
      vtkIdType *cPtr(nullptr), *dPtr(nullptr);
      vtkSmartPointer<vtkIdTypeArray> cellLocations(vtkSmartPointer<vtkIdTypeArray>::New());
      {
        cellLocations->SetNumberOfComponents(1);
        cellLocations->SetNumberOfTuples(nbCells);
        cPtr = cellLocations->GetPointer(0);
      }
      vtkSmartPointer<vtkIdTypeArray> cells(vtkSmartPointer<vtkIdTypeArray>::New());
      {
        cells->SetNumberOfComponents(1);
        cells->SetNumberOfTuples(mVor->getNodalConnectivity()->getNumberOfTuples());
        dPtr = cells->GetPointer(0);
      }
      const mcIdType *connPtr(mVor->getNodalConnectivity()->begin()),
        *connIPtr(mVor->getNodalConnectivityIndex()->begin());
      for (int i = 0; i < nbCells; i++, connIPtr++)
      {
        *dPtr++ = 2;
        dPtr = std::copy(connPtr + connIPtr[0] + 1, connPtr + connIPtr[1], dPtr);
        *cPtr++ = 3 * i;
      }
      vtkSmartPointer<vtkCellArray> cells2(vtkSmartPointer<vtkCellArray>::New());
      cells2->SetCells(nbCells, cells);
      ret->SetCells(cellTypes, cellLocations, cells2);
      break;
    }
    default:
      throw INTERP_KERNEL::Exception("Not implemented yet !");
  }
  ret->SetPoints(points);
  return ret;
}

MCAuto<DataArrayDouble> ForceBuilder(const std::vector<std::size_t>& TAB, const DataArrayDouble* matrix, const DataArrayDouble* eqn)
{
  MCAuto<DataArrayDouble> tmp0, tmp1, ret;
  tmp0 = matrix->keepSelectedComponents({ TAB[0] });
  tmp1 = eqn->keepSelectedComponents({ 0 });
  MCAuto<DataArrayDouble> p0(DataArrayDouble::Multiply(tmp0, tmp1));
  tmp0 = matrix->keepSelectedComponents({ TAB[1] });
  tmp1 = eqn->keepSelectedComponents({ 1 });
  MCAuto<DataArrayDouble> p1(DataArrayDouble::Multiply(tmp0, tmp1));
  ret = DataArrayDouble::Add(p0, p1);
  tmp0 = matrix->keepSelectedComponents({ TAB[2] });
  tmp1 = eqn->keepSelectedComponents({ 2 });
  MCAuto<DataArrayDouble> p2(DataArrayDouble::Multiply(tmp0, tmp1));
  ret = DataArrayDouble::Add(ret, p2);
  return ret;
}

double ReturnInertia(
  const double pOut[3], const DataArrayDouble* OM, const DataArrayDouble* area_field_ids)
{
  MCAuto<DataArrayDouble> base_X(DataArrayDouble::New());
  base_X->alloc(OM->getNumberOfTuples(), 3);
  base_X->setPartOfValuesSimple1(pOut[0], 0, OM->getNumberOfTuples(), 1, 0, 1, 1);
  base_X->setPartOfValuesSimple1(pOut[1], 0, OM->getNumberOfTuples(), 1, 1, 2, 1);
  base_X->setPartOfValuesSimple1(pOut[2], 0, OM->getNumberOfTuples(), 1, 2, 3, 1);
  MCAuto<DataArrayDouble> dist_to_base_X;
  {
    MCAuto<DataArrayDouble> tmp(DataArrayDouble::Dot(OM, base_X));
    tmp = DataArrayDouble::Multiply(tmp, base_X);
    tmp = DataArrayDouble::Substract(OM, tmp);
    dist_to_base_X = tmp->magnitude();
  }
  MCAuto<DataArrayDouble> inertiaArr;
  {
    MCAuto<DataArrayDouble> tmp(DataArrayDouble::Multiply(dist_to_base_X, dist_to_base_X));
    inertiaArr = DataArrayDouble::Multiply(tmp, area_field_ids);
  }
  double inertiaTmp;
  inertiaArr->accumulate(&inertiaTmp);
  return inertiaTmp;
}

void FindPrincipalAxeInternal(const double startVector[3], const double normalFace[3],
  const DataArrayDouble* OM, const DataArrayDouble* area_field_ids,
  const std::vector<double>& posToIterate, double& angleDegree, double outputAxis[3],
  double& inertia)
{
  constexpr double CENTER[3] = { 0, 0, 0 };
  inertia = -std::numeric_limits<double>::max();
  for (auto zePos : posToIterate)
  {
    double p[3] = { startVector[0], startVector[1], startVector[2] }, pOut[3];
    DataArrayDouble::Rotate3DAlg(CENTER, normalFace, zePos / 180. * M_PI, 1, p, pOut);
    double inertiaTmp = ReturnInertia(pOut, OM, area_field_ids);
    if (inertiaTmp > inertia)
    {
      inertia = inertiaTmp;
      std::copy(pOut, pOut + 3, outputAxis);
      angleDegree = zePos;
    }
  }
}

void FindPrincipalAxe(const double startVector[3], const double normalFace[3],
  const DataArrayDouble* OM, const DataArrayDouble* area_field_ids, double& angleDegree,
  double outputAxis[3], double& inertia)
{
  std::vector<double> R(
    { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170 });
  FindPrincipalAxeInternal(
    startVector, normalFace, OM, area_field_ids, R, angleDegree, outputAxis, inertia);
  for (int i = 0; i < 5; ++i)
  {
    std::vector<double> Q({ -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 });
    const double CST(std::pow((double)10., (double)-i));
    std::for_each(Q.begin(), Q.end(), [angleDegree, CST](double& v) { v = CST * v + angleDegree; });
    FindPrincipalAxeInternal(
      startVector, normalFace, OM, area_field_ids, Q, angleDegree, outputAxis, inertia);
  }
}

vtkSmartPointer<vtkTable> ComputeTorseurCIH(vtkUnstructuredGrid* usgIn)
{
  std::vector<MCAuto<MEDCouplingUMesh> > m;
  {
    std::vector<MCAuto<DataArrayIdType> > ids;
    ConvertFromUnstructuredGrid(usgIn, m, ids);
  }
  vtkDataArray* sief(nullptr);
  {
    int nArrays(usgIn->GetPointData()->GetNumberOfArrays());
    for (int i = 0; i < nArrays; i++)
    {
      vtkDataArray* array(usgIn->GetPointData()->GetArray(i));
      if (!array)
        continue;
      std::string name(array->GetName());
      if (name.find("SIEF") != std::string::npos)
        if (array->GetNumberOfComponents() == 6)
        {
          if (sief)
          {
            std::ostringstream oss;
            oss << "ComputeTorseurCIH : several candidates for SIEF field !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
          sief = array;
        }
    }
  }
  if (!sief)
    throw INTERP_KERNEL::Exception("ComputeTorseurCIH : unable to find a field for SIEF!");
  MCAuto<MEDCouplingFieldDouble> area_field(m[0]->getMeasureField(true));
  double area;
  area_field->accumulate(&area); // 1
  MCAuto<DataArrayDouble> centerOfMassField(m[0]->computeCellCenterOfMass());
  double centerOfMass[3];
  {
    MCAuto<DataArrayDouble> tmp(
      DataArrayDouble::Multiply(centerOfMassField, area_field->getArray()));
    tmp->accumulate(centerOfMass);
    std::for_each(centerOfMass, centerOfMass + 3, [area](double& v) { v /= area; });
  } // 2
  m[0]->unPolyze();
  std::set<INTERP_KERNEL::NormalizedCellType> s(m[0]->getAllGeoTypes());
  if (s.size() != 1)
  {
    std::ostringstream oss;
    oss << "ComputeTorseurCIH : Only TRI3 are supported !";
    throw INTERP_KERNEL::Exception(oss.str());
  }
  if (*(s.begin()) != INTERP_KERNEL::NORM_TRI3)
  {
    std::ostringstream oss;
    oss << "ComputeTorseurCIH : Only TRI3 are supported !";
    throw INTERP_KERNEL::Exception(oss.str());
  }
  MCAuto<MEDCouplingFieldDouble> f(MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES));
  {
    f->setMesh(m[0]);
    MCAuto<DataArrayDouble> tmp(ConvertVTKArrayToMCArrayDouble(sief));
    f->setArray(tmp);
  }
  MCAuto<MEDCouplingFieldDouble> fCell(f->nodeToCellDiscretization());
  MCAuto<DataArrayIdType> ids;
  {
    MCAuto<DataArrayIdType> tmp(area_field->getArray()->findIdsLowerThan(1e-7));
    ids = tmp->buildComplement(m[0]->getNumberOfCells());
  }
  MCAuto<MEDCouplingUMesh> m_ids(m[0]->buildPartOfMySelf(ids->begin(), ids->end()));
  MCAuto<DataArrayDouble> eqn;
  {
    MCAuto<DataArrayDouble> tmp(m_ids->computePlaneEquationOf3DFaces());
    eqn = tmp->keepSelectedComponents({ 0, 1, 2 });
    tmp = eqn->magnitude();
    eqn = DataArrayDouble::Divide(eqn, tmp);
  }
  MCAuto<DataArrayDouble> area_field_ids(
    area_field->getArray()->selectByTupleId(ids->begin(), ids->end()));
  MCAuto<DataArrayDouble> area_vector(DataArrayDouble::Multiply(eqn, area_field_ids));
  MCAuto<DataArrayDouble> matrix(fCell->getArray()->selectByTupleId(ids->begin(), ids->end()));
  MCAuto<DataArrayDouble> F_x, F_y, F_z;
  {
    F_x = ForceBuilder({ 0, 3, 4 }, matrix, eqn);
    F_y = ForceBuilder({ 3, 1, 5 }, matrix, eqn);
    F_z = ForceBuilder({ 4, 5, 2 }, matrix, eqn);
  }
  //
  //
  MCAuto<DataArrayDouble> F(DataArrayDouble::Meld({ F_x, F_y, F_z }));
  F->multiplyEqual(area_vector);
  double ZeForce[3], normalFace[3];
  F->accumulate(ZeForce);
  eqn->accumulate(normalFace);
  {
    double normalFaceNorm(sqrt(normalFace[0] * normalFace[0] + normalFace[1] * normalFace[1] +
      normalFace[2] * normalFace[2]));
    std::for_each(normalFace, normalFace + 3, [normalFaceNorm](double& v) { v /= normalFaceNorm; });
  }
  double ForceNormale[3]; // 3
  {
    double tmp(
      ZeForce[0] * normalFace[0] + ZeForce[1] * normalFace[1] + ZeForce[2] * normalFace[2]);
    ForceNormale[0] = tmp * normalFace[0];
    ForceNormale[1] = tmp * normalFace[1];
    ForceNormale[2] = tmp * normalFace[2];
  }
  double TangentForce[3] = { ZeForce[0] - ForceNormale[0], ZeForce[1] - ForceNormale[1],
    ZeForce[2] - ForceNormale[2] }; // 4
  MCAuto<DataArrayDouble> bary(m_ids->computeCellCenterOfMass());
  MCAuto<DataArrayDouble> OM;
  {
    MCAuto<DataArrayDouble> centerOfMass2(DataArrayDouble::New());
    centerOfMass2->alloc(1, 3);
    std::copy(centerOfMass, centerOfMass + 3, centerOfMass2->getPointer());
    OM = DataArrayDouble::Substract(bary, centerOfMass2);
  }
  double momentum[3]; // 5
  {
    MCAuto<DataArrayDouble> tmp(DataArrayDouble::CrossProduct(OM, F));
    tmp->accumulate(momentum);
  }
  double InertiaNormale(ReturnInertia(normalFace, OM, area_field_ids)); // 6
  double base[9];
  DataArrayDouble::GiveBaseForPlane(normalFace, base);
  double angleDegree, outputAxis[3], inertia;
  FindPrincipalAxe(base, normalFace, OM, area_field_ids, angleDegree, outputAxis, inertia);
  double tangentOther[3] = { normalFace[1] * outputAxis[2] - normalFace[2] * outputAxis[1],
    normalFace[2] * outputAxis[0] - normalFace[0] * outputAxis[2],
    normalFace[0] * outputAxis[1] - normalFace[1] * outputAxis[0] };
  double inertiaOther(ReturnInertia(tangentOther, OM, area_field_ids));
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
  // std::cerr << "########################################## vtkTorseurCIH::RequestInformation
  // ##########################################" << std::endl;
  try
  {
    vtkSmartPointer<vtkUnstructuredGrid> usgIn;
    ExtractInfo(inputVector[0], usgIn);
  }
  catch (INTERP_KERNEL::Exception& e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkTorseurCIH::RequestInformation : " << e.what()
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

int vtkTorseurCIH::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // std::cerr << "########################################## vtkTorseurCIH::RequestData
  // ##########################################" << std::endl;
  try
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
  catch (INTERP_KERNEL::Exception& e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkTorseurCIH::RequestData : " << e.what() << std::endl;
    if (this->HasObserver("ErrorEvent"))
      this->InvokeEvent("ErrorEvent", const_cast<char*>(oss.str().c_str()));
    else
      vtkOutputWindowDisplayErrorText(const_cast<char*>(oss.str().c_str()));
    vtkObject::BreakOnError();
    return 0;
  }
  return 1;
}

void vtkTorseurCIH::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
