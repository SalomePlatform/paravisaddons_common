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

#ifndef vtkRateOfFlowThroughSection_h__
#define vtkRateOfFlowThroughSection_h__

#include <vtkDataObjectAlgorithm.h>
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
#include <vtkSetGet.h>

#include "VTKToMEDMem.h"

#include <map>
#include <deque>
#include <sstream>

class vtkMutableDirectedGraph;

class VTK_EXPORT vtkRateOfFlowThroughSection : public vtkDataObjectAlgorithm
{
public:
  static vtkRateOfFlowThroughSection *New();
  vtkTypeMacro(vtkRateOfFlowThroughSection, vtkDataObjectAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  void SetSourceData(vtkDataObject *input);

  void SetSourceConnection(vtkAlgorithmOutput *algOutput);

  int FillOutputPortInformation(int, vtkInformation *) override;

protected:
  vtkRateOfFlowThroughSection();
  ~vtkRateOfFlowThroughSection() override;

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  int NumberOfTimeSteps;
  int CurrentTimeIndex;
  bool IsExecuting;
  class  VTK_EXPORT vtkInternal {
  public:
    vtkInternal() {};
    void pushData(double timeStep, double value) { _data.emplace_back(timeStep, value); }
    void fillTable(vtkTable *table) const;
    void analyzeInputDataSets(vtkUnstructuredGrid *ds1, vtkDataSet *ds2);
    bool computationNeeded() const
    {
      if (_recomputationOfMatrixNeeded)
      {
        _matrix.clear();
      }
      return _recomputationOfMatrixNeeded;
    };
    std::vector<std::map<int, double>> &getMatrix() { return _matrix; }
    void setOrtho(const MEDCoupling::DataArrayDouble *ortho) { _ortho.takeRef(ortho); }
    const MEDCoupling::DataArrayDouble *getOrtho() const { return (const MEDCoupling::DataArrayDouble *)_ortho; }
    void setMeasure(const MEDCoupling::DataArrayDouble *measure) { _measure.takeRef(measure); }
    const MEDCoupling::DataArrayDouble *getMeasure() { return (const MEDCoupling::DataArrayDouble *)_measure; }

  private:
    std::vector<std::pair<double, double>> _data;
    vtkMTimeType _mt1 = 0;
    vtkMTimeType _mt2 = 0;
    mutable std::vector<std::map<int, double>> _matrix;
    bool _recomputationOfMatrixNeeded = true;
    MEDCoupling::MCConstAuto<MEDCoupling::DataArrayDouble> _ortho;
    MEDCoupling::MCConstAuto<MEDCoupling::DataArrayDouble> _measure;
  };

  vtkInternal *Internal;

private:
  vtkRateOfFlowThroughSection(const vtkRateOfFlowThroughSection &) = delete;
  void operator=(const vtkRateOfFlowThroughSection &) = delete;
};

#endif
