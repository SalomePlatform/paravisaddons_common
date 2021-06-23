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

#include "vtkContactReader.h"

#include <vtkArrowSource.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDelimitedTextReader.h>
#include <vtkDoubleArray.h>
#include <vtkErrorCode.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPVGlyphFilter.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTable.h>
#include <vtkVariant.h>
#include <vtkVariantArray.h>

#include <exception>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>

constexpr char INST_TAG[] = "INST";

class MyException : public std::exception
{
public:
  MyException(const char *what) : _what(what) {}
  MyException(const std::string &what) : _what(what) {}
  ~MyException() throw() {}
  const char *what() const throw() { return _what.c_str(); }

private:
  std::string _what;
};

template <class T>
class AutoPtr
{
public:
  AutoPtr(T *ptr = 0) : _ptr(ptr) {}
  ~AutoPtr() { destroyPtr(); }
  bool isNull() const { return _ptr == 0; }
  bool isNotNull() const { return !isNull(); }
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

bool PresenceOf(const std::vector<std::string> &arr, const std::string &what)
{
  auto pos = std::find(arr.begin(), arr.end(), what);
  return pos!=arr.end();
}

vtkIdType PosOf(const std::vector<std::string> &arr, const std::string &what)
{
  auto pos = std::find(arr.begin(), arr.end(), what);
  if (pos == arr.end())
  {
    std::ostringstream oss;
    oss << "vtkContactReader::PosOf : Fail to locate \"" << what << "\" in array !";
    throw MyException(oss.str());
  }
  return std::distance(arr.begin(), pos);
}

double InstValueOf(vtkTable *table, vtkIdType iRow, vtkIdType iCol)
{
  vtkVariantArray *row(table->GetRow(iRow));
  vtkVariant *elt(row->GetPointer(iCol));
  return elt->ToDouble();
}

vtkSmartPointer<vtkTable> LoadDataFromFile(const char *FileName)
{
  vtkNew<vtkDelimitedTextReader> reader;
  reader->SetFileName(FileName);
  reader->SetDetectNumericColumns(true);
  reader->SetUseStringDelimiter(true);
  reader->SetHaveHeaders(true);
  reader->SetFieldDelimiterCharacters(" ");
  reader->SetAddTabFieldDelimiter(true);
  reader->SetMergeConsecutiveDelimiters(true);
  reader->Update();
  vtkTable *table(reader->GetOutput());
  vtkSmartPointer<vtkTable> ret(table);
  return ret;
}

std::vector<std::string> GetColumnNames(vtkTable *table)
{
  vtkIdType nbCols(table->GetNumberOfColumns());
  std::vector<std::string> colNames(nbCols);
  for (vtkIdType iCol = 0; iCol < nbCols; iCol++)
  {
    colNames[iCol] = table->GetColumnName(iCol);
  }
  return colNames;
}

void FillValue(vtkVariantArray *row, double *ptToFeed, std::size_t ipos, vtkIdType pos)
{
  bool isOK(false);
  vtkVariant *elt(row->GetPointer(pos));
  ptToFeed[ipos] = elt->ToDouble(&isOK);
  if (!isOK)
  {
    std::ostringstream oss;
    oss << "vtkContactReader::FillValue : Error during analyze content of file ! Float64 expected !";
    throw MyException(oss.str());
  }
}

vtkStandardNewMacro(vtkContactReader);

vtkContactReader::vtkContactReader() : FileName(NULL), ScaleFactor(0.02),IsTimed(false)
{
  this->SetNumberOfInputPorts(0);
}

vtkContactReader::~vtkContactReader()
{
}

int vtkContactReader::RequestInformation(vtkInformation *vtkNotUsed(request),
                                         vtkInformationVector **vtkNotUsed(inputVector),
                                         vtkInformationVector *outputVector)
{
  vtkInformation *outInfo(outputVector->GetInformationObject(0));
  //
  vtkSmartPointer<vtkTable> table(LoadDataFromFile(this->FileName));
  std::vector<std::string> colNames(GetColumnNames(table));
  this->IsTimed = PresenceOf(colNames,INST_TAG);
  vtkIdType nbRows(table->GetNumberOfRows());
  if(!this->IsTimed)
    return 1;
  vtkIdType InstPos = PosOf(colNames, INST_TAG);
  std::set<double> allInsts;
  for (vtkIdType iRow = 0; iRow < nbRows; iRow++)
  {
    vtkVariantArray *row(table->GetRow(iRow));
    vtkVariant *elt(row->GetPointer(InstPos));
    allInsts.insert(elt->ToDouble());
  }
  std::vector<double> allInstsV(allInsts.begin(),allInsts.end());
  double timeRange[2];
  timeRange[0]=allInstsV.front();
  timeRange[1]=allInstsV.back();
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),allInstsV.data(),(int)allInstsV.size());
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),timeRange,2);
  return 1;
}


int vtkContactReader::RequestData(vtkInformation *vtkNotUsed(request),
                                  vtkInformationVector **vtkNotUsed(inputVector),
                                  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo(outputVector->GetInformationObject(0));
  vtkPolyData *output(vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
  //
  try
  {
    vtkSmartPointer<vtkTable> table(LoadDataFromFile(this->FileName));
    vtkIdType nbRows(table->GetNumberOfRows()), nbCols(table->GetNumberOfColumns());
    std::vector<std::string> colNames(GetColumnNames(table));
    vtkIdType XPos(PosOf(colNames, "X")), YPos(PosOf(colNames, "Y")), ZPos(PosOf(colNames, "Z")), DXPos(PosOf(colNames, "DX")), DYPos(PosOf(colNames, "DY")), DZPos(PosOf(colNames, "DZ"));
    // check of presence of INST. If yes -> it means that input file is on different time steps.
    vtkIdType InstPos(std::numeric_limits<vtkIdType>::max());
    //
    vtkSmartPointer<vtkDoubleArray> coords(vtkSmartPointer<vtkDoubleArray>::New()), vectArr(vtkSmartPointer<vtkDoubleArray>::New());
    //
    vtkIdType nbTuples(nbRows),startRow(0);
    double reqTS(0.0);
    if(this->IsTimed)
    {
      InstPos = PosOf(colNames, INST_TAG);
      if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
        reqTS=outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
      vtkIdType iRow = 0;
      for (; iRow < nbRows && InstValueOf(table,iRow,InstPos)!=reqTS; iRow++);
      startRow = iRow;
      for (; iRow < nbRows && InstValueOf(table,iRow,InstPos)==reqTS; iRow++);
      nbTuples = iRow-startRow;
    }
    //
    vectArr->SetNumberOfComponents(3);
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(nbTuples);
    vectArr->SetNumberOfTuples(nbTuples);
    double *ptToFeed1(coords->Begin()), *ptToFeed2(vectArr->Begin());
    const vtkIdType POS[3] = {XPos, YPos, ZPos}, DX[3] = {DXPos, DYPos, DZPos};
    for (vtkIdType iRow = startRow; iRow < startRow+nbTuples; iRow++, ptToFeed1 += 3, ptToFeed2 += 3)
    {
      vtkVariantArray *row(table->GetRow(iRow));
      for (std::size_t ipos = 0; ipos < 3; ipos++)
      {
        FillValue(row, ptToFeed1, ipos, POS[ipos]);
        FillValue(row, ptToFeed2, ipos, DX[ipos]);
      }
      std::for_each(ptToFeed2, ptToFeed2 + 3, [](double &v) { v = -v; });
    }
    vectArr->SetName("Resultante");
    vtkNew<vtkPolyData> ret;
    vtkSmartPointer<vtkPoints> pts(vtkSmartPointer<vtkPoints>::New());
    pts->SetData(coords);
    ret->SetPoints(pts);
    ret->GetPointData()->AddArray(vectArr);
    //
    vtkNew<vtkPVGlyphFilter> glyph;
    glyph->SetInputData(ret);
    glyph->SetGlyphMode(0);       //vtkPVGlyphFilter::ALL_POINTS
    glyph->SetVectorScaleMode(0); //vtkPVGlyphFilter::SCALE_BY_MAGNITUDE
    //
    vtkNew<vtkArrowSource> arrow;
    arrow->SetTipResolution(6);
    arrow->SetTipRadius(0.1);
    arrow->SetTipLength(0.35);
    arrow->SetShaftResolution(6);
    arrow->SetShaftRadius(0.03);
    glyph->SetSourceConnection(arrow->GetOutputPort());
    //idx,port,connection,fieldAssociation,name
    glyph->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Resultante"); //idx==0 -> scaleArray
    glyph->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Resultante"); //idx==1 -> orientationArray
    glyph->SetScaleFactor(this->ScaleFactor);
    glyph->Update();
    output->ShallowCopy(glyph->GetOutput());
    if(this->IsTimed)
      output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),reqTS);
    output->GetPointData()->SetActiveAttribute(0, vtkDataSetAttributes::SCALARS);
    //output->ShallowCopy(ret);
  }
  catch (MyException &e)
  {
    vtkErrorMacro(<< "vtkContactReader::RequestData : during read of " << this->FileName << " : " << e.what());
    return 0;
  }
  return 1;
}

void vtkContactReader::SetScaleFactor(double newScaleFactor)
{
  if (this->ScaleFactor != newScaleFactor)
  {
    this->ScaleFactor = newScaleFactor;
    this->Modified();
  }
}

void vtkContactReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
