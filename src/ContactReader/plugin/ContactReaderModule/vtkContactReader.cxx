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

vtkStandardNewMacro(vtkContactReader);

vtkContactReader::vtkContactReader() : FileName(NULL), ScaleFactor(0.02)
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
  return 1;
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

int vtkContactReader::RequestData(vtkInformation *vtkNotUsed(request),
                                  vtkInformationVector **vtkNotUsed(inputVector),
                                  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo(outputVector->GetInformationObject(0));
  vtkPolyData *output(vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
  //
  try
  {
    vtkNew<vtkDelimitedTextReader> reader;
    reader->SetFileName(this->FileName);
    reader->SetDetectNumericColumns(true);
    reader->SetUseStringDelimiter(true);
    reader->SetHaveHeaders(true);
    reader->SetFieldDelimiterCharacters(" ");
    reader->SetAddTabFieldDelimiter(true);
    reader->SetMergeConsecutiveDelimiters(true);
    reader->Update();
    vtkTable *table(reader->GetOutput());
    vtkIdType nbRows(table->GetNumberOfRows()), nbCols(table->GetNumberOfColumns());
    std::vector<std::string> colNames(nbCols);
    for (vtkIdType iCol = 0; iCol < nbCols; iCol++)
    {
      colNames[iCol] = table->GetColumnName(iCol);
    }
    vtkIdType XPos(PosOf(colNames, "X")), YPos(PosOf(colNames, "Y")), ZPos(PosOf(colNames, "Z")), DXPos(PosOf(colNames, "DX")), DYPos(PosOf(colNames, "DY")), DZPos(PosOf(colNames, "DZ"));
    //
    vtkSmartPointer<vtkDoubleArray> coords(vtkSmartPointer<vtkDoubleArray>::New()), vectArr(vtkSmartPointer<vtkDoubleArray>::New());
    vectArr->SetNumberOfComponents(3);
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(nbRows);
    vectArr->SetNumberOfTuples(nbRows);
    double *ptToFeed1(coords->Begin()), *ptToFeed2(vectArr->Begin());
    const vtkIdType POS[3] = {XPos, YPos, ZPos}, DX[3] = {DXPos, DYPos, DZPos};
    for (vtkIdType iRow = 0; iRow < nbRows; iRow++, ptToFeed1 += 3, ptToFeed2 += 3)
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
