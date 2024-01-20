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
// Author : Anthony Geay (EDF R&D)

#ifndef vtkSedimentDeposit_h__
#define vtkSedimentDeposit_h__

#include <vtkDataObjectAlgorithm.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSetGet.h>
#include <vtkDataSet.h>
#include <vtkTable.h>
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

#include <vector>
#include <deque>
#include <map>
#include <sstream>
#include <memory>

class VTK_EXPORT vtkSedimentDeposit : public vtkDataObjectAlgorithm
{
public:
  static vtkSedimentDeposit *New();
  vtkTypeMacro(vtkSedimentDeposit, vtkDataObjectAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  void SetSourceData(vtkDataObject *input);
  void SetSourceConnection(vtkAlgorithmOutput *algOutput);

  int FillOutputPortInformation(int, vtkInformation *) override;

protected:
  vtkSedimentDeposit();
  ~vtkSedimentDeposit() override = default;

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  int NumberOfTimeSteps;
  int CurrentTimeIndex;
  bool IsExecuting;
  class  VTK_EXPORT vtkInternal
  {
  public:
    vtkInternal(int curveId, int nbOfCurves) :_curve_id(curveId), _nb_of_curves(nbOfCurves) {}
    void pushData(double timeStep, double positiveValue, double negativeValue) { _data.emplace_back(timeStep, positiveValue, negativeValue); }
    void fillTable(vtkTable *table) const;
    void analyzeInputDataSets(vtkUnstructuredGrid *ds1, vtkDataSet *ds2);
    bool computationNeeded() const;
    MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> &meshOrigin() { return _meshorig; }
    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> &untouched2DCells() { return _untouched_2d_cells; }
    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> &cellsAtBoundary() { return _cells_at_boundary_origin; }
    MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> &centers() { return _centers; }
    MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> &measure() { return _vol; }
    std::string getReprDependingPos(const std::string& origName) const;

  private:
    std::vector<std::tuple<double, double, double>> _data;
    vtkMTimeType _mt1 = 0;
    vtkMTimeType _mt2 = 0;
    int _curve_id;
    int _nb_of_curves;
    bool _recomputationOfMatrixNeeded = true;
    mutable MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> _meshorig;
    mutable MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> _untouched_2d_cells;
    mutable MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> _cells_at_boundary_origin;
    mutable MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> _centers;
    mutable MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> _vol;
  };
  std::vector< std::unique_ptr< vtkInternal > > Internal2;

private:
  vtkSedimentDeposit(const vtkSedimentDeposit &) = delete;
  void operator=(const vtkSedimentDeposit &) = delete;
};

#endif
