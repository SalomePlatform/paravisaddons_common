// Copyright (C) 2021-2023  CEA, EDF
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

#include "vtkQuadraticToLinear.h"

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
#include <vtkLongArray.h>
#include <vtkExecutive.h>
#include <vtkFloatArray.h>
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
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkTimeStamp.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVariantArray.h>
#include <vtkWarpScalar.h>

#include <deque>
#include <map>
#include <sstream>

vtkStandardNewMacro(vtkQuadraticToLinear);

constexpr int NB_QUAD_TO_LINEAR_1 = 7;

constexpr VTKCellType QUAD_TO_LINEAR_QUAD1[NB_QUAD_TO_LINEAR_1] = {VTK_QUADRATIC_EDGE, VTK_QUADRATIC_TRIANGLE, VTK_QUADRATIC_QUAD, VTK_QUADRATIC_TETRA, VTK_QUADRATIC_WEDGE, VTK_QUADRATIC_PYRAMID, VTK_QUADRATIC_HEXAHEDRON};

constexpr VTKCellType QUAD_TO_LINEAR_LIN1[NB_QUAD_TO_LINEAR_1] = {VTK_LINE, VTK_TRIANGLE, VTK_QUAD, VTK_TETRA, VTK_WEDGE, VTK_PYRAMID, VTK_HEXAHEDRON};

constexpr int NB_QUAD_TO_LINEAR_2 = 4;

constexpr VTKCellType QUAD_TO_LINEAR_QUAD2[NB_QUAD_TO_LINEAR_2] = {VTK_BIQUADRATIC_TRIANGLE, VTK_BIQUADRATIC_QUAD, VTK_BIQUADRATIC_QUADRATIC_WEDGE, VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON};

constexpr VTKCellType QUAD_TO_LINEAR_LIN2[NB_QUAD_TO_LINEAR_2] = {VTK_TRIANGLE, VTK_QUAD, VTK_WEDGE, VTK_HEXAHEDRON};

constexpr int NB_NB_NODES_PER_CELL = 7;

constexpr VTKCellType NB_NODES_PER_CELL_1[NB_NB_NODES_PER_CELL] = {VTK_LINE, VTK_TRIANGLE, VTK_QUAD, VTK_TETRA, VTK_WEDGE, VTK_PYRAMID, VTK_HEXAHEDRON};

constexpr int NB_NODES_PER_CELL_2[NB_NB_NODES_PER_CELL] = {2, 3, 4, 4, 6, 5, 8};

///////////////////

template <class T>
class AutoPtr
{
public:
  AutoPtr(T *ptr = 0) : _ptr(ptr) {}
  ~AutoPtr() { destroyPtr(); }
  AutoPtr &operator=(T *ptr)
  {
    if (_ptr != ptr)
    {
      destroyPtr();
      _ptr = ptr;
    }
    return *this;
  }
  T *operator->() { return _ptr; }
  const T *operator->() const { return _ptr; }
  T &operator*() { return *_ptr; }
  const T &operator*() const { return *_ptr; }
  operator T *() { return _ptr; }
  operator const T *() const { return _ptr; }

private:
  void destroyPtr() { delete[] _ptr; }

private:
  T *_ptr;
};

class MZCException : public std::exception
{
public:
  MZCException(const std::string &s) : _reason(s) {}
  virtual const char *what() const throw() { return _reason.c_str(); }
  virtual ~MZCException() throw() {}

private:
  std::string _reason;
};

void ExtractInfo(vtkInformationVector *inputVector, vtkUnstructuredGrid *&usgIn)
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
      throw MZCException("Input dataSet must be a DataSet or single elt multi block dataset expected !");
    if (input1->GetNumberOfBlocks() != 1)
      throw MZCException("Input dataSet is a multiblock dataset with not exactly one block ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    vtkDataObject *input2(input1->GetBlock(0));
    if (!input2)
      throw MZCException("Input dataSet is a multiblock dataset with exactly one block but this single element is NULL !");
    vtkDataSet *input2c(vtkDataSet::SafeDownCast(input2));
    if (!input2c)
      throw MZCException("Input dataSet is a multiblock dataset with exactly one block but this single element is not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    input = input2c;
  }
  if (!input)
    throw MZCException("Input data set is NULL !");
  usgIn = vtkUnstructuredGrid::SafeDownCast(input);
  if (!usgIn)
    throw MZCException("Input data set is not an unstructured mesh ! This filter works only on unstructured meshes !");
}

vtkSmartPointer<vtkDataArray> Reduce(const int *new2Old, int newNbPts, vtkDataArray *array)
{
  if (!array)
    throw MZCException("Reduce : null input vector !");
  int nbOfCompo(array->GetNumberOfComponents());
  vtkSmartPointer<vtkDataArray> zeRet;
  if (vtkDoubleArray::SafeDownCast(array))
  {
    vtkSmartPointer<vtkFloatArray> ret(vtkSmartPointer<vtkFloatArray>::New());
    zeRet = ret;
    ret->SetNumberOfComponents(nbOfCompo);
    ret->SetNumberOfTuples(newNbPts);
    vtkDoubleArray *array1(vtkDoubleArray::SafeDownCast(array));
    if (array1)
    {
      const double *inpCoords(array1->GetPointer(0));
      float *outCoords(ret->GetPointer(0));
      for (int i = 0; i < newNbPts; i++, outCoords += nbOfCompo)
        std::copy(inpCoords + new2Old[i] * nbOfCompo, inpCoords + (new2Old[i] + 1) * nbOfCompo, outCoords);
    }
    else
    {
      std::ostringstream oss;
      oss << "Only Double array managed for the moment in input !" << array->GetName();
      throw MZCException(oss.str());
    }
  }
  else if (vtkFloatArray::SafeDownCast(array))
  {
    vtkSmartPointer<vtkFloatArray> ret(vtkSmartPointer<vtkFloatArray>::New());
    zeRet = ret;
    ret->SetNumberOfComponents(nbOfCompo);
    ret->SetNumberOfTuples(newNbPts);
    vtkFloatArray *array1(vtkFloatArray::SafeDownCast(array));
    if (array1)
    {
      const float *inpCoords(array1->GetPointer(0));
      float *outCoords(ret->GetPointer(0));
      for (int i = 0; i < newNbPts; i++, outCoords += nbOfCompo)
        std::copy(inpCoords + new2Old[i] * nbOfCompo, inpCoords + (new2Old[i] + 1) * nbOfCompo, outCoords);
    }
    else
    {
      std::ostringstream oss;
      oss << "Only Float array managed for the moment in input !" << array->GetName();
      throw MZCException(oss.str());
    }
  }
  else if (vtkLongArray::SafeDownCast(array))
  {
    vtkSmartPointer<vtkLongArray> ret(vtkSmartPointer<vtkLongArray>::New());
    zeRet = ret;
    ret->SetNumberOfComponents(nbOfCompo);
    ret->SetNumberOfTuples(newNbPts);
    vtkLongArray *array1(vtkLongArray::SafeDownCast(array));
    if (array1)
    {
      const long *inpCoords(array1->GetPointer(0));
      long *outCoords(ret->GetPointer(0));
      for (int i = 0; i < newNbPts; i++, outCoords += nbOfCompo)
        std::copy(inpCoords + new2Old[i] * nbOfCompo, inpCoords + (new2Old[i] + 1) * nbOfCompo, outCoords);
    }
    else
    {
      std::ostringstream oss;
      oss << "Only Long array managed for the moment in input !" << array->GetName();
      throw MZCException(oss.str());
    }
  }
  else if (vtkIntArray::SafeDownCast(array))
  {
    vtkSmartPointer<vtkIntArray> ret(vtkSmartPointer<vtkIntArray>::New());
    zeRet = ret;
    ret->SetNumberOfComponents(nbOfCompo);
    ret->SetNumberOfTuples(newNbPts);
    vtkIntArray *array1(vtkIntArray::SafeDownCast(array));
    if (array1)
    {
      const int *inpCoords(array1->GetPointer(0));
      int *outCoords(ret->GetPointer(0));
      for (int i = 0; i < newNbPts; i++, outCoords += nbOfCompo)
        std::copy(inpCoords + new2Old[i] * nbOfCompo, inpCoords + (new2Old[i] + 1) * nbOfCompo, outCoords);
    }
    else
    {
      std::ostringstream oss;
      oss << "Only int32 array managed for the moment in input !" << array->GetName();
      throw MZCException(oss.str());
    }
  }
  else
  {
    std::ostringstream oss;
    oss << "Reduce : unmanaged type ! Check field " << array->GetName();
    throw MZCException(oss.str());
  }
  for (int i = 0; i < nbOfCompo; i++)
  {
    const char *compoName(array->GetComponentName(i));
    zeRet->SetComponentName(i, compoName);
  }
  return zeRet;
}

vtkSmartPointer<vtkUnstructuredGrid> ComputeQuadToLinear(vtkUnstructuredGrid *usg)
{
  std::map<VTKCellType, VTKCellType> linToQuad;
  std::map<VTKCellType, vtkIdType> nbNodesPerType;
  for (int i = 0; i < NB_QUAD_TO_LINEAR_1; i++)
    linToQuad[QUAD_TO_LINEAR_QUAD1[i]] = QUAD_TO_LINEAR_LIN1[i];
  for (int i = 0; i < NB_QUAD_TO_LINEAR_2; i++)
    linToQuad[QUAD_TO_LINEAR_QUAD2[i]] = QUAD_TO_LINEAR_LIN2[i];
  for (int i = 0; i < NB_NB_NODES_PER_CELL; i++)
    nbNodesPerType[NB_NODES_PER_CELL_1[i]] = NB_NODES_PER_CELL_2[i];
  if (!usg)
    throw MZCException("ComputeQuadToLinear : null input pointer !");
  vtkIdType nbPts = usg->GetNumberOfPoints();
  vtkIdType nbOfCells = usg->GetNumberOfCells();
  vtkIdType maxNbOfNodesPerCell = 0;
  std::vector<bool> old2NewVB(nbPts, false);
  for (vtkIdType i = 0; i < nbOfCells; i++)
  {
    vtkCell *cell(usg->GetCell(i));
    VTKCellType ct((VTKCellType)cell->GetCellType());
    std::map<VTKCellType, VTKCellType>::const_iterator it(linToQuad.find(ct));
    if (it != linToQuad.end())
    {
      std::map<VTKCellType, vtkIdType>::const_iterator it2(nbNodesPerType.find((*it).second));
      maxNbOfNodesPerCell = std::max(it2->second, maxNbOfNodesPerCell);
      for (vtkIdType j = 0; j < (*it2).second; j++)
      {
        vtkIdType ptId(cell->GetPointId(j));
        old2NewVB[ptId] = true;
      }
    }
    else
    {
      if (ct != VTK_POLYHEDRON)
      {
        vtkIdType nbPtsOfCell = cell->GetNumberOfPoints();
        maxNbOfNodesPerCell = std::max(nbPtsOfCell, maxNbOfNodesPerCell);
        for (vtkIdType j = 0; j < nbPtsOfCell; j++)
        {
          vtkIdType ptId = cell->GetPointId(j);
          old2NewVB[ptId] = true;
        }
      }
      else
        throw MZCException("ComputeQuadToLinear : polyhedrons are not managed yet !");
    }
  }
  int newNbPts(std::count(old2NewVB.begin(), old2NewVB.end(), true));
  AutoPtr<int> new2Old(new int[newNbPts]), old2New(new int[nbPts]);
  struct Renumberer
  {
    Renumberer(int *pt) : _cnt(0), _pt(pt) {}
    void operator()(bool v)
    {
      if (v)
      {
        *(_pt++) = _cnt;
      }
      _cnt++;
    }

  private:
    int _cnt;
    int *_pt;
  };
  struct RenumbererR
  {
    RenumbererR(int *pt) : _cnt(0), _pt(pt) {}
    void operator()(bool v)
    {
      *_pt++ = _cnt;
      if (v)
        _cnt++;
    }

  private:
    int _cnt;
    int *_pt;
  };
  std::for_each(old2NewVB.begin(), old2NewVB.end(), Renumberer(new2Old));
  std::for_each(old2NewVB.begin(), old2NewVB.end(), RenumbererR(old2New));
  //
  vtkNew<vtkUnstructuredGrid> ret;
  vtkNew<vtkPoints> pts;
  // deal with coordinates
  vtkSmartPointer<vtkDataArray> newCoords(Reduce(new2Old, newNbPts, usg->GetPoints()->GetData()));
  //
  ret->Initialize();
  ret->Allocate(nbOfCells);
  { // deal with connectivity
    AutoPtr<vtkIdType> connOfCellTmp(new vtkIdType[maxNbOfNodesPerCell]);
    for (vtkIdType i = 0; i < nbOfCells; i++)
    {
      vtkCell *cell(usg->GetCell(i));
      VTKCellType ct((VTKCellType)cell->GetCellType());
      std::map<VTKCellType, VTKCellType>::const_iterator it(linToQuad.find(ct));
      if (it != linToQuad.end())
      {
        std::map<VTKCellType, vtkIdType>::const_iterator it2(nbNodesPerType.find(it->second));
        for (vtkIdType j = 0; j < it2->second; j++)
          connOfCellTmp[j] = old2New[cell->GetPointId(j)];
        ret->InsertNextCell(it->second, it2->second, connOfCellTmp);
      }
      else
      {
        vtkIdType nbPtsOfCell(cell->GetNumberOfPoints());
        for (vtkIdType j = 0; j < nbPtsOfCell; j++)
          connOfCellTmp[j] = old2New[cell->GetPointId(j)];
        ret->InsertNextCell(ct, nbPtsOfCell, connOfCellTmp);
      }
    }
  }
  ret->SetPoints(pts);
  pts->SetData(newCoords);
  // Deal with cell fields
  if (usg->GetCellData())
  {
    ret->GetCellData()->ShallowCopy(usg->GetCellData());
  }
  if (usg->GetPointData())
  {
    vtkPointData *pd(usg->GetPointData());
    for (int i = 0; i < pd->GetNumberOfArrays(); i++)
    {
      vtkSmartPointer<vtkDataArray> arr(Reduce(new2Old, newNbPts, pd->GetArray(i)));
      arr->SetName(pd->GetArray(i)->GetName());
      ret->GetPointData()->AddArray(arr);
    }
  }
  //
  return ret;
}

////////////////////

int vtkQuadraticToLinear::RequestInformation(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkQuadraticToLinear::RequestInformation ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid *usgIn(0);
    ExtractInfo(inputVector[0], usgIn);
  }
  catch (MZCException &e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkQuadraticToLinear::RequestInformation : " << e.what() << std::endl;
    if (this->HasObserver("ErrorEvent"))
      this->InvokeEvent("ErrorEvent", const_cast<char *>(oss.str().c_str()));
    else
      vtkOutputWindowDisplayErrorText(const_cast<char *>(oss.str().c_str()));
    vtkObject::BreakOnError();
    return 0;
  }
  return 1;
}

int vtkQuadraticToLinear::RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  //std::cerr << "########################################## vtkQuadraticToLinear::RequestData        ##########################################" << std::endl;
  try
  {
    vtkUnstructuredGrid *usgIn(0);
    ExtractInfo(inputVector[0], usgIn);
    vtkSmartPointer<vtkUnstructuredGrid> ret(ComputeQuadToLinear(usgIn));
    vtkInformation *outInfo(outputVector->GetInformationObject(0));
    vtkUnstructuredGrid *output(vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
    output->ShallowCopy(ret);
  }
  catch (MZCException &e)
  {
    std::ostringstream oss;
    oss << "Exception has been thrown in vtkQuadraticToLinear::RequestInformation : " << e.what() << std::endl;
    if (this->HasObserver("ErrorEvent"))
      this->InvokeEvent("ErrorEvent", const_cast<char *>(oss.str().c_str()));
    else
      vtkOutputWindowDisplayErrorText(const_cast<char *>(oss.str().c_str()));
    vtkObject::BreakOnError();
    return 0;
  }
  return 1;
}

void vtkQuadraticToLinear::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
/*
(cd /home/H87074/salome/DEV/tools/build/Paravisaddons-master-cm302-pv501 ; export CURRENT_SOFTWARE_SRC_DIR=/home/H87074/salome/DEV/tools/src/PARAVISADDONS ; export CURRENT_SOFTWARE_BUILD_DIR=/home/H87074/salome/DEV/tools/build/Paravisaddons-master-cm302-pv501 ; export CURRENT_SOFTWARE_INSTALL_DIR=/home/H87074/salome/DEV/tools/install/Paravisaddons-master-cm302-pv501 ; export PYTHON_VERSION="2.7" ; . /home/H87074/salome/DEV/salome_modules.sh >/dev/null 2>&1 ; . /home/H87074/salome/DEV/tools/build/Paravisaddons-master-cm302-pv501/.yamm/env_build.sh >/dev/null 2>&1 ; . /home/H87074/salome/DEV/salome_prerequisites.sh >/dev/null 2>&1 ; cmake -DCMAKE_INSTALL_PREFIX=/home/H87074/salome/DEV/tools/install/Paravisaddons-master-cm302-pv501   -DCMAKE_MODULE_PATH=${CONFIGURATION_CMAKE_DIR}   -DCMAKE_BUILD_TYPE=Debug  /home/H87074/salome/DEV/tools/src/PARAVISADDONS )

(cd /home/H87074/salome/DEV/tools/build/Paravisaddons-master-cm302-pv501/QuadraticToLinear ; export CURRENT_SOFTWARE_SRC_DIR=/home/H87074/salome/DEV/tools/src/PARAVISADDONS ; export CURRENT_SOFTWARE_BUILD_DIR=/home/H87074/salome/DEV/tools/build/Paravisaddons-master-cm302-pv501 ; export CURRENT_SOFTWARE_INSTALL_DIR=/home/H87074/salome/DEV/tools/install/Paravisaddons-master-cm302-pv501 ; export PYTHON_VERSION="2.7" ; . /home/H87074/salome/DEV/salome_modules.sh >/dev/null 2>&1 ; . /home/H87074/salome/DEV/tools/build/Paravisaddons-master-cm302-pv501/.yamm/env_build.sh >/dev/null 2>&1 ; . /home/H87074/salome/DEV/salome_prerequisites.sh >/dev/null 2>&1 ; make -j8 install )

/home/H87074/salome/prerequisites/src/ParaView/VTK/Common/DataModel/vtkCellType.h:87
*/
