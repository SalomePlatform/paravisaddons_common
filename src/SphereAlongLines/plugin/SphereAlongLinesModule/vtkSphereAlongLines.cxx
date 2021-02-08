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

#include "vtkSphereAlongLines.h"

#include <vtkAlgorithmOutput.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkExecutive.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <deque>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>

vtkStandardNewMacro(vtkSphereAlongLines);

template <class T>
struct VTKTraits
{
};

template <>
struct VTKTraits<float>
{
  using VtkType = vtkFloatArray;
};

template <>
struct VTKTraits<double>
{
  using VtkType = vtkDoubleArray;
};

constexpr const char INTEGRATIONTIME_ARR_NAME[] = "IntegrationTime";

///////////////////

class MZCException : public std::exception
{
public:
  MZCException(const std::string& s)
    : _reason(s)
  {
  }
  virtual const char* what() const throw() { return _reason.c_str(); }
  virtual ~MZCException() throw() {}

private:
  std::string _reason;
};

void ExtractInfo(vtkInformationVector* inputVector, vtkPolyData*& usgIn)
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
      throw MZCException(
        "Input dataSet must be a DataSet or single elt multi block dataset expected !");
    if (input1->GetNumberOfBlocks() != 1)
      throw MZCException("Input dataSet is a multiblock dataset with not exactly one block ! Use "
                         "MergeBlocks or ExtractBlocks filter before calling this filter !");
    vtkDataObject* input2(input1->GetBlock(0));
    if (!input2)
      throw MZCException("Input dataSet is a multiblock dataset with exactly one block but this "
                         "single element is NULL !");
    vtkDataSet* input2c(vtkDataSet::SafeDownCast(input2));
    if (!input2c)
      throw MZCException(
        "Input dataSet is a multiblock dataset with exactly one block but this single element is "
        "not a dataset ! Use MergeBlocks or ExtractBlocks filter before calling this filter !");
    input = input2c;
  }
  if (!input)
    throw MZCException("Input data set is NULL !");
  vtkPointData* att(input->GetPointData());
  vtkPolyData* zeInput(vtkPolyData::SafeDownCast(input));
  if (!zeInput)
    throw MZCException("Input dataSet is not a polydata as expected !");
  usgIn = zeInput;
}

vtkDataArray* GetCoords(vtkPointSet* ds)
{
  vtkPoints* pts(ds->GetPoints());
  if (!pts)
    throw MZCException("GetCoords : internal error 2 !");
  vtkDataArray* data(pts->GetData());
  if (!data)
    throw MZCException("GetCoords : internal error 3 !");
  if (data->GetNumberOfComponents() != 3)
    throw MZCException("GetCoords : internal error 4 !");
  return data;
}

vtkDoubleArray* GetIntegrationTime(vtkPointSet* ds)
{
  vtkDataSetAttributes* dsa(ds->GetPointData());
  if (!dsa)
    throw MZCException(
      "GetIntegrationTime : no point data ! Is the input data comes from Stream Tracer !");
  int idx(0);
  vtkAbstractArray* arr(dsa->GetAbstractArray(INTEGRATIONTIME_ARR_NAME, idx));
  if (!arr)
  {
    std::ostringstream oss;
    oss << "GetIntegrationTime : no such " << INTEGRATIONTIME_ARR_NAME
        << " array in input dataset !";
    throw MZCException(oss.str());
  }
  vtkDoubleArray* ret(vtkDoubleArray::SafeDownCast(arr));
  if (!ret)
  {
    std::ostringstream oss;
    oss << "GetIntegrationTime :" << INTEGRATIONTIME_ARR_NAME << " array expected to be float64 !";
    throw MZCException(oss.str());
  }
  if (ret->GetNumberOfComponents() != 1)
  {
    std::ostringstream oss;
    oss << "GetIntegrationTime :" << INTEGRATIONTIME_ARR_NAME
        << " array expected to be single compo !";
    throw MZCException(oss.str());
  }
  return ret;
}

template <class T>
void RearrangeIfNecessaryImpl(typename VTKTraits<T>::VtkType* coords, vtkDoubleArray* intTime,
  std::vector<std::vector<vtkIdType> >& connect, double eps)
{
  std::size_t nbCells(connect.size());
  if (nbCells % 2 != 0)
    return;
  std::size_t nbCellsCand(nbCells / 2);
  const double* intTimePtr(intTime->GetPointer(0));
  const T* coordsPtr(coords->GetPointer(0));
  T distance[3];
  std::vector<std::vector<vtkIdType> > connectOut(nbCellsCand);
  for (std::size_t i = 0; i < nbCellsCand; i++)
  {
    const std::vector<vtkIdType>& elt0(connect[i]);
    const std::vector<vtkIdType>& elt1(connect[i + nbCellsCand]);
    if (elt0.empty() || elt1.empty())
      return;
    vtkIdType pt0(elt0[0]), pt1(elt1[0]);
    if (pt0 == pt1)
      continue;
    if (std::abs(intTimePtr[pt0]) > eps || std::abs(intTimePtr[pt1]) > eps)
      return;
    std::transform(coordsPtr + 3 * pt0, coordsPtr + 3 * (pt0 + 1), coordsPtr + 3 * pt1, distance,
      [](T a, T b) { return b - a; });
    if (std::sqrt(
          distance[0] * distance[0] + distance[1] * distance[1] + distance[2] * distance[2]) > eps)
      return;
    std::vector<vtkIdType>& tab(connectOut[i]);
    tab.insert(tab.end(), elt1.rbegin(), elt1.rend());
    tab.back() = elt0[0];
    tab.insert(tab.end(), elt0.begin() + 1, elt0.end());
  }
  connect = connectOut;
}

void RearrangeIfNecessary(vtkDataArray* coords, vtkDoubleArray* intTime,
  std::vector<std::vector<vtkIdType> >& connect, double eps)
{
  vtkFloatArray* c0(vtkFloatArray::SafeDownCast(coords));
  if (c0)
  {
    RearrangeIfNecessaryImpl<float>(c0, intTime, connect, eps);
    return;
  }
  vtkDoubleArray* c1(vtkDoubleArray::SafeDownCast(coords));
  if (c1)
  {
    RearrangeIfNecessaryImpl<double>(c1, intTime, connect, eps);
    return;
  }
  throw MZCException("Not recognized type of data for coordinates !");
}

class vtkSphereAlongLines::vtkInternals
{
  friend class CurvAbsPathWalker;

  class PathWalker
  {
  protected:
    PathWalker(const vtkSphereAlongLines::vtkInternals* intern)
      : _internal(intern)
    {
    }

  public:
    virtual double getTimeFrom(std::size_t pathId, double absTime2) const = 0;
    virtual bool isInside(std::size_t pathId, double realTime, vtkIdType& a, vtkIdType& b,
      double& aw, double& bw) const = 0;

  protected:
    const vtkSphereAlongLines::vtkInternals* _internal;
  };

  class IntTimeGlobPathWalker : public PathWalker
  {
  public:
    IntTimeGlobPathWalker(const vtkSphereAlongLines::vtkInternals* intern)
      : PathWalker(intern)
    {
    }
    double getTimeFrom(std::size_t pathId, double absTime2) const;
    bool isInside(std::size_t pathId, double realTime, vtkIdType& a, vtkIdType& b, double& aw,
      double& bw) const override;
  };

  class LocalPathWalker : public PathWalker
  {
  protected:
    LocalPathWalker(const vtkSphereAlongLines::vtkInternals* intern)
      : PathWalker(intern)
    {
    }
    bool isInside(std::size_t pathId, double realTime, vtkIdType& a, vtkIdType& b, double& aw,
      double& bw) const override;
    double getTimeFrom(std::size_t pathId, double absTime2) const override;
    virtual const std::vector<double>& getRankingArray(std::size_t pathId) const = 0;
  };

  class IntTimeLocPathWalker : public LocalPathWalker
  {
  public:
    IntTimeLocPathWalker(const vtkSphereAlongLines::vtkInternals* intern)
      : LocalPathWalker(intern)
    {
    }
    const std::vector<double>& getRankingArray(std::size_t pathId) const override
    {
      return _internal->_integration_time_along_pathes[pathId];
    }
  };

  class CurvAbsPathWalker : public LocalPathWalker
  {
  public:
    CurvAbsPathWalker(const vtkSphereAlongLines::vtkInternals* intern)
      : LocalPathWalker(intern)
    {
    }
    const std::vector<double>& getRankingArray(std::size_t pathId) const override
    {
      return _internal->_curv_absc_along_pathes[pathId];
    }
  };

public:
  void initCacheIfNeeded(vtkPolyData* ds);
  vtkSmartPointer<vtkPolyData> dataSetAtNormalizedTime(
    vtkPolyData* ds, double absTime, int walkType) const;
  std::unique_ptr<PathWalker> buildPathWalker(int walkType) const;

private:
  void initCacheForce(vtkPolyData* ds);

private:
  vtkMTimeType _input_DS_time = 0;
  std::vector<std::vector<vtkIdType> > _connectivity;
  //! for each path it stores the integration time the corresponding time. Every array is expected
  //! to be in ascending order
  std::vector<std::vector<double> > _integration_time_along_pathes;
  //! for each path it stores curv abscissa along pathes
  std::vector<std::vector<double> > _curv_absc_along_pathes;
  std::vector<std::pair<double, double> > _time_range_per_path;
  std::vector<std::pair<double, double> > _ca_range_per_path;
  double _abs_min = std::numeric_limits<double>::max();
  double _abs_max = -std::numeric_limits<double>::max();
};

double vtkSphereAlongLines::vtkInternals::IntTimeGlobPathWalker::getTimeFrom(
  std::size_t, double absTime2) const
{
  return _internal->_abs_min + absTime2 * (_internal->_abs_max - _internal->_abs_min);
}

bool vtkSphereAlongLines::vtkInternals::IntTimeGlobPathWalker::isInside(
  std::size_t pathId, double realTime, vtkIdType& a, vtkIdType& b, double& aw, double& bw) const
{
  const std::pair<double, double>& pathRange(_internal->_time_range_per_path[pathId]);
  if (realTime < pathRange.first || realTime > pathRange.second)
    return false;
  const std::vector<double>& timeAlongPath(_internal->_integration_time_along_pathes[pathId]);
  const std::vector<vtkIdType>& associatedNodes(_internal->_connectivity[pathId]);
  std::vector<double>::const_iterator it(std::find_if(timeAlongPath.begin(), timeAlongPath.end(),
    [realTime](const double& val) -> bool { return realTime < val; }));
  std::size_t pos(std::distance(timeAlongPath.begin(), it));
  pos = std::min(pos, timeAlongPath.size() - 1);
  if (*it != realTime)
  {
    a = associatedNodes[pos - 1];
    b = associatedNodes[pos];
    aw = (timeAlongPath[pos] - realTime) / (timeAlongPath[pos] - timeAlongPath[pos - 1]);
    bw = 1. - aw;
  }
  else
  {
    a = associatedNodes[pos];
    b = a;
    aw = 1.;
    bw = 0.;
  }
  return true;
}

bool vtkSphereAlongLines::vtkInternals::LocalPathWalker::isInside(
  std::size_t pathId, double realTime, vtkIdType& a, vtkIdType& b, double& aw, double& bw) const
{
  const std::vector<double>& timeAlongPath(getRankingArray(pathId));
  const std::vector<vtkIdType>& associatedNodes(_internal->_connectivity[pathId]);
  std::vector<double>::const_iterator it(std::find_if(timeAlongPath.begin(), timeAlongPath.end(),
    [realTime](const double& val) -> bool { return realTime <= val; }));
  std::size_t pos(std::distance(timeAlongPath.begin(), it));
  pos = std::min(pos, timeAlongPath.size() - 1);
  if (*it != realTime)
  {
    a = associatedNodes[pos - 1];
    b = associatedNodes[pos];
    aw = (timeAlongPath[pos] - realTime) / (timeAlongPath[pos] - timeAlongPath[pos - 1]);
    bw = 1. - aw;
  }
  else
  {
    a = associatedNodes[pos];
    b = a;
    aw = 1.;
    bw = 0.;
  }
  return true;
}

double vtkSphereAlongLines::vtkInternals::LocalPathWalker::getTimeFrom(
  std::size_t pathId, double absTime2) const
{
  double realTime(std::numeric_limits<double>::max());
  {
    const std::vector<double>& timeAlongPath(getRankingArray(pathId));
    realTime = timeAlongPath.front() + absTime2 * (timeAlongPath.back() - timeAlongPath.front());
    realTime = std::max(realTime, timeAlongPath.front());
    realTime = std::min(realTime, timeAlongPath.back());
  }
  return realTime;
}

std::unique_ptr<vtkSphereAlongLines::vtkInternals::PathWalker>
vtkSphereAlongLines::vtkInternals::buildPathWalker(int walkType) const
{
  switch (walkType)
  {
    case 0:
      return std::unique_ptr<vtkSphereAlongLines::vtkInternals::PathWalker>(
        new IntTimeGlobPathWalker(this));
    case 1:
      return std::unique_ptr<vtkSphereAlongLines::vtkInternals::PathWalker>(
        new IntTimeLocPathWalker(this));
    case 2:
      return std::unique_ptr<vtkSphereAlongLines::vtkInternals::PathWalker>(
        new CurvAbsPathWalker(this));
    default:
      return std::unique_ptr<vtkSphereAlongLines::vtkInternals::PathWalker>(
        new CurvAbsPathWalker(this));
  }
}

void vtkSphereAlongLines::vtkInternals::initCacheIfNeeded(vtkPolyData* ds)
{
  vtkMTimeType mtime(ds->GetMTime());
  if (ds->GetMTime() == _input_DS_time)
    return;
  initCacheForce(ds);
  _input_DS_time = mtime;
}

void vtkSphereAlongLines::vtkInternals::initCacheForce(vtkPolyData* ds)
{
  // std::cout << "Force cache" << std::endl;
  _abs_min = std::numeric_limits<double>::max();
  _abs_max = -std::numeric_limits<double>::max();
  _connectivity.clear();
  _integration_time_along_pathes.clear();
  _curv_absc_along_pathes.clear();
  _time_range_per_path.clear();
  _ca_range_per_path.clear();
  // to Improve
  vtkDataArray* cooInBase(GetCoords(ds));
  vtkFloatArray* cooIn(vtkFloatArray::SafeDownCast(cooInBase));
  //
  if ((ds->GetPolys() && ds->GetPolys()->GetNumberOfCells() > 0) ||
    (ds->GetStrips() && ds->GetStrips()->GetNumberOfCells() > 0) ||
    (ds->GetVerts() && ds->GetVerts()->GetNumberOfCells() > 0))
    throw MZCException(
      "Presence of strips/vertices/polygons in input polydata ! Invalid with this type of filter!");
  vtkCellArray* cc(ds->GetLines());
  if (!cc)
    throw MZCException("No polylines in input polydata ! Difficult to build something on it !");
  vtkIdType nbCells(cc->GetNumberOfCells());
  _connectivity.resize(nbCells);
  vtkIdType npts;
  const vtkIdType* pts;
  cc->InitTraversal();
  for (vtkIdType i = 0; cc->GetNextCell(npts, pts); i++)
  {
    std::vector<vtkIdType>& conn2(_connectivity[i]);
    conn2.insert(conn2.end(), pts, pts + npts);
  }
  cc->InitTraversal();
  vtkDoubleArray* ita(GetIntegrationTime(ds));
  RearrangeIfNecessary(GetCoords(ds), ita, _connectivity, 1e-5);
  const double* itaPtr(ita->GetPointer(0));
  std::size_t nbCellsEff(_connectivity.size());
  _integration_time_along_pathes.resize(nbCellsEff);
  _curv_absc_along_pathes.resize(nbCellsEff);
  _time_range_per_path.resize(nbCellsEff);
  _ca_range_per_path.resize(nbCellsEff);
  for (std::size_t i = 0; i < nbCellsEff; i++)
  {
    const std::vector<vtkIdType>& conn2(_connectivity[i]);
    std::vector<double>& vals(_integration_time_along_pathes[i]);
    std::vector<double>& absCurv(_curv_absc_along_pathes[i]);
    std::transform(conn2.begin(), conn2.end(), std::back_inserter(vals),
      [itaPtr](vtkIdType nodeId) { return itaPtr[nodeId]; });
    float lastVal(0.);
    const float* cooInPtr(cooIn->GetPointer(conn2[0] * 3));
    std::transform(conn2.begin(), conn2.end(), std::back_inserter(absCurv),
      [cooIn, &cooInPtr, &lastVal](vtkIdType nodeId) {
        const float* pt(cooIn->GetPointer(nodeId * 3));
        lastVal += sqrt((cooInPtr[0] - pt[0]) * (cooInPtr[0] - pt[0]) +
          (cooInPtr[1] - pt[1]) * (cooInPtr[1] - pt[1]) +
          (cooInPtr[2] - pt[2]) * (cooInPtr[2] - pt[2]));
        cooInPtr = pt;
        return lastVal;
      });
    _time_range_per_path[i].first = vals.front();
    _time_range_per_path[i].second = vals.back();
    _ca_range_per_path[i].first = absCurv.front();
    _ca_range_per_path[i].second = absCurv.back();
    _abs_min = std::min(_abs_min, vals.front());
    _abs_max = std::max(_abs_max, vals.back());
    // std::cout << "cell " << i << " -> "; std::for_each(conn2.begin(),conn2.end(),[](vtkIdType
    // elt) { std::cout << elt << " "; }); std::cout << std::endl; std::cout << "cell " << i << " ->
    // "; std::for_each(absCurv.begin(),absCurv.end(),[](double elt) { std::cout << elt << " "; });
    // std::cout << std::endl;
  }
  // std::cerr << _abs_min << " - " << _abs_max << std::endl;
}

template <class VTKDATAARRAY>
void dealWith(
  VTKDATAARRAY* arr1, vtkDataArray* arrOut, vtkIdType a, double aw, vtkIdType b, double bw)
{
  int nbCompo(arr1->GetNumberOfComponents());
  using VT = typename VTKDATAARRAY::ValueType;
  VT *tupleA(arr1->GetPointer(nbCompo * a)), *tupleB(arr1->GetPointer(b * nbCompo));
  {
    VT* tmpData = new VT[nbCompo];
    std::transform(tupleA, tupleA + nbCompo, tupleB, tmpData,
      [aw, bw](VT a, VT b) -> VT { return VT(a) * VT(aw) + VT(b) * VT(bw); });
    arrOut->InsertNextTuple(tmpData);
    delete[] tmpData;
  }
}

class CooAssign
{
public:
  virtual void apply(
    vtkIdType a, double aw, vtkIdType b, double bw, vtkDoubleArray* cooOut) const = 0;
};

template <class T>
class CooAssignT : public CooAssign
{
public:
  using VtkType = typename VTKTraits<T>::VtkType;
  CooAssignT(VtkType* arr)
    : _coo_in(arr->GetPointer(0))
  {
  }
  void apply(vtkIdType a, double aw, vtkIdType b, double bw, vtkDoubleArray* cooOut) const override
  {
    const T *pt0(_coo_in + 3 * a), *pt1(_coo_in + 3 * b);
    std::transform(pt0, pt0 + 3, pt1, ptToAdd, [aw, bw](T elt0, T elt1) -> double {
      return double(elt0) * double(aw) + double(elt1) * double(bw);
    });
    cooOut->InsertNextTuple(ptToAdd);
  }

private:
  const T* _coo_in;
  T tmp[3];
  mutable double ptToAdd[3];
};

vtkSmartPointer<vtkPolyData> vtkSphereAlongLines::vtkInternals::dataSetAtNormalizedTime(
  vtkPolyData* ds, double absTime, int walkType) const
{
  double absTime2(std::min(std::max(absTime, 0.), 1.));
  std::size_t maxNbOfPts(_integration_time_along_pathes.size());
  std::vector<vtkSmartPointer<vtkDataArray> > outArrays;
  std::vector<vtkDataArray*> inArrays;
  {
    vtkDataSetAttributes* dsa(ds->GetPointData());
    for (int i = 0; i < dsa->GetNumberOfArrays(); i++)
    {
      vtkDataArray* arr(dsa->GetArray(i));
      if (!arr)
        continue;
      vtkDoubleArray* arr1(vtkDoubleArray::SafeDownCast(arr));
      if (arr1)
      {
        vtkSmartPointer<vtkDoubleArray> outArray(vtkDoubleArray::New());
        outArray->SetName(arr->GetName());
        outArray->SetNumberOfComponents(arr->GetNumberOfComponents());
        outArray->SetNumberOfTuples(0);
        outArrays.push_back(outArray);
        inArrays.push_back(arr);
        continue;
      }
      vtkFloatArray* arr2(vtkFloatArray::SafeDownCast(arr));
      if (arr2)
      {
        vtkSmartPointer<vtkFloatArray> outArray(vtkFloatArray::New());
        outArray->SetName(arr->GetName());
        outArray->SetNumberOfComponents(arr->GetNumberOfComponents());
        outArray->SetNumberOfTuples(0);
        outArrays.push_back(outArray);
        inArrays.push_back(arr);
        continue;
      }
    }
  }
  std::size_t nbArrays(outArrays.size());
  vtkDataArray* cooInBase(GetCoords(ds));
  std::unique_ptr<CooAssign> cooInDS;
  vtkFloatArray* cooIn(vtkFloatArray::SafeDownCast(cooInBase));
  if (cooIn)
  {
    std::unique_ptr<CooAssignT<float> > tmp(new CooAssignT<float>(cooIn));
    cooInDS = std::move(tmp);
  }
  vtkDoubleArray* cooIn2(vtkDoubleArray::SafeDownCast(cooInBase));
  if (cooIn2)
  {
    std::unique_ptr<CooAssignT<double> > tmp(new CooAssignT<double>(cooIn2));
    cooInDS = std::move(tmp);
  }
  const float* cooInPtr(cooIn->GetPointer(0));
  // TODO : improve
  vtkSmartPointer<vtkPolyData> outDS(vtkPolyData::New());
  vtkSmartPointer<vtkDoubleArray> cooOut(vtkDoubleArray::New());
  cooOut->SetNumberOfComponents(3);
  cooOut->SetNumberOfTuples(0);
  {
    vtkNew<vtkPoints> pts;
    pts->SetData(cooOut);
    outDS->SetPoints(pts);
  }
  vtkNew<vtkCellArray> verts;
  vtkNew<vtkIdTypeArray> conn;
  conn->SetNumberOfComponents(1);
  vtkIdType nbOfPts(0);
  vtkIdType connToAdd[2];
  connToAdd[0] = 1;
  std::vector<vtkIdType> connv;
  std::unique_ptr<vtkSphereAlongLines::vtkInternals::PathWalker> pw(buildPathWalker(walkType));
  for (std::size_t i = 0; i < maxNbOfPts; i++)
  {
    vtkIdType a, b;
    double aw, bw;
    double realTime(pw->getTimeFrom(i, absTime2));
    if (pw->isInside(i, realTime, a, b, aw, bw))
    {
      // std::cout << i << " - " << a << " - " << b << " * " << aw << " * " << bw << std::endl;
      cooInDS->apply(a, aw, b, bw, cooOut);
      connToAdd[1] = nbOfPts++;
      connv.insert(connv.end(), connToAdd, connToAdd + 2);
      for (std::size_t j = 0; j < nbArrays; j++)
      {
        vtkDataArray* arr(inArrays[j]);
        vtkDoubleArray* arr1(vtkDoubleArray::SafeDownCast(arr));
        if (arr1)
        {
          dealWith(arr1, outArrays[j], a, aw, b, bw);
          continue;
        }
        vtkFloatArray* arr2(vtkFloatArray::SafeDownCast(arr));
        if (arr2)
        {
          dealWith(arr2, outArrays[j], a, aw, b, bw);
          continue;
        }
      }
    }
  }
  conn->SetNumberOfTuples(connv.size());
  std::copy(connv.begin(), connv.end(), conn->GetPointer(0));
  verts->SetCells(nbOfPts, conn);
  outDS->SetVerts(verts);
  {
    vtkDataSetAttributes* dsa(outDS->GetPointData());
    for (auto arr : outArrays)
    {
      dsa->AddArray(arr);
    }
  }
  return outDS;
}

////////////////////

vtkSphereAlongLines::vtkSphereAlongLines()
  : Internal(new vtkInternals)
  , AnimationTime(0.)
  , WalkType(2)
{
}

vtkSphereAlongLines::~vtkSphereAlongLines()
{
  delete this->Internal;
}

int vtkSphereAlongLines::RequestInformation(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // std::cerr << "##########################################
  // vtkSphereAlongLines::RequestInformation ##########################################" <<
  // std::endl;
  try
  {
    vtkPolyData* ds(nullptr);
    ExtractInfo(inputVector[0], ds);
  }
  catch (MZCException& e)
  {
    vtkErrorMacro(<< "Exception has been thrown in vtkSphereAlongLines::RequestInformation : "
                  << e.what());
    return 0;
  }
  return 1;
}

int vtkSphereAlongLines::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // std::cerr << "########################################## vtkSphereAlongLines::RequestData
  // ##########################################" << std::endl;
  try
  {
    vtkPolyData* ds(nullptr);
    ExtractInfo(inputVector[0], ds);
    this->Internal->initCacheIfNeeded(ds);
    vtkSmartPointer<vtkPolyData> pts(
      this->Internal->dataSetAtNormalizedTime(ds, this->AnimationTime, this->WalkType));
    vtkInformation* outInfo(outputVector->GetInformationObject(0));
    vtkPolyData* output(vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
    output->ShallowCopy(pts);
  }
  catch (MZCException& e)
  {
    vtkErrorMacro(<< "Exception has been thrown in vtkSphereAlongLines::Data : " << e.what());
    return 0;
  }
  return 1;
}

void vtkSphereAlongLines::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
