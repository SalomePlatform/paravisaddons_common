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

#include "vtkExplodePolyLine.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkMergeBlocks.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLineSource.h>
#include <vtkMultiBlockDataGroupFilter.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkRibbonFilter.h>
#include <vtkPVGlyphFilter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkTable.h>
#include <vtkTessellatorFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVariant.h>
#include <vtkVariantArray.h>
#include <vtkMeshQuality.h>
#include <vtkCellCenters.h>

#include <cmath>
#include <sstream>
vtkStandardNewMacro(vtkExplodePolyLine);

class MyException : public std::exception
{
public:
  MyException(const std::string& s):_reason(s) { }
  virtual const char *what() const throw() { return _reason.c_str(); }
  virtual ~MyException() throw() { }
private:
  std::string _reason;
};

//-----------------------------------------------------------------------------
void vtkExplodePolyLine::ExtractInfo(vtkInformationVector* inputVector, vtkSmartPointer<vtkPolyData>& usgIn)
{
  vtkInformation* inputInfo(inputVector->GetInformationObject(0));
  vtkPolyData* input(vtkPolyData::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
  if (!input)
  {
    vtkUnstructuredGrid* input2(vtkUnstructuredGrid::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT())));
    if(!input2)
      vtkErrorMacro("Input data set is not vtkPolyData !");

    vtkNew<vtkDataSetSurfaceFilter> surface;

    surface->SetNonlinearSubdivisionLevel(0);
    surface->SetInputData(input2);
    surface->Update();
    usgIn = surface->GetOutput();
    return;
  }
  usgIn = vtkPolyData::SafeDownCast(input);
}


int vtkExplodePolyLine::FillInputPortInformation(int vtkNotUsed(port), vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

int vtkExplodePolyLine::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int vtkExplodePolyLine::RequestData(vtkInformation* vtkNotUsed(request),vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  try
  {
    vtkInformation* outInfo(outputVector->GetInformationObject(0));
    vtkPointSet* output(vtkPointSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
    vtkSmartPointer<vtkPolyData> usgIn;
    this->ExtractInfo(inputVector[0], usgIn);
    vtkCellArray *cc(usgIn->GetLines());
    vtkIdType nbEntries(cc->GetNumberOfConnectivityEntries());
    const vtkIdType *conn(cc->GetData()->GetPointer(0));
    vtkNew<vtkPolyData> pd;
    vtkNew<vtkCellArray> cb;
    if(nbEntries<1)
    {
      throw MyException("Input PolyData looks empty !");
    }
    if(conn[0] == nbEntries-1)
    {
      vtkIdType nbCells(nbEntries-2);
      vtkNew<vtkIdTypeArray> connOut;
      connOut->SetNumberOfComponents(1);
      connOut->SetNumberOfTuples(3*nbCells);
      vtkIdType *connOutPtr(connOut->GetPointer(0));
      for( auto iCell = 0 ; iCell < nbCells ; ++iCell, connOutPtr+=3 )
      {
        connOutPtr[0] = 2;
        connOutPtr[1] = conn[1+iCell];
        connOutPtr[2] = conn[1+iCell+1];
      }
      cb->SetCells(nbCells,connOut);
    }
    else if(nbEntries == 3*cc->GetNumberOfCells())
    {
      output->ShallowCopy(usgIn);
      return 1;
    }
    else
    {
      throw MyException("Input PolyData is not containing only lines as required as precondition !");
    }
    pd->SetLines(cb);
    vtkNew<vtkPoints> pts;
    // here DeepCopy is required, because GetMeshMTime of input PolyData is modified and generate useless computation afterwards in the pipeline. Bug ?
    pts->DeepCopy(usgIn->GetPoints());
    pd->SetPoints(pts);
    output->ShallowCopy(pd);
  }
  catch(MyException& e)
  {
    vtkErrorMacro("Exception has been thrown in vtkComplexMode::RequestInformation : " << e.what());
  }
  return 1;
}

//-----------------------------------------------------------------------------
void vtkExplodePolyLine::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
