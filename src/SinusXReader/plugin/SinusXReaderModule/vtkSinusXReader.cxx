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

#include "vtkSinusXReader.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkErrorCode.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkStreamingDemandDrivenPipeline.h>

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
  AutoPtr(T *ptr = nullptr) : _ptr(ptr) {}
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

bool isFloat(const std::string &st, double &val)
{
  if (st == "NaN")
  {
    val = 0.;
    return true;
  }
  std::istringstream iss(st);
  iss >> val;
  return iss.eof() && !iss.fail() && !iss.bad();
}

bool FindPart(const std::string &st, std::size_t &work, double &val)
{
  std::string part0;
  std::size_t pos0(st.find_first_not_of(" \t", work));
  if (pos0 == std::string::npos)
    return false;
  std::size_t pos1(st.find_first_of(" \t", pos0));
  if (pos1 == std::string::npos)
    return false;
  part0 = st.substr(pos0, pos1 - pos0);
  if (!isFloat(part0, val))
    return false;
  work = pos1;
  return true;
}

bool isDataLine(const char *line, std::streamsize szOfLine, double &val0, double &val1, double &val2)
{
  if (szOfLine < 1)
    return false;
  std::string st(line, szOfLine - 1);
  std::size_t work(0);
  if (!FindPart(st, work, val0))
    return false;
  if (!FindPart(st, work, val1))
    return false;
  if (!FindPart(st, work, val2))
    return false;
  std::size_t pos0(st.find_first_not_of(" \t", work));
  if (pos0 == std::string::npos)
    return false;
  std::size_t pos1(st.find_first_of(" \t", pos0));
  if (pos1 != std::string::npos)
    return false;
  std::string endOfLine(st.substr(pos0));
  if (endOfLine.length() != 1)
    return false;
  const char c(endOfLine[0]);
  if (c < 'A' || c > 'Z')
    return false;
  return true;
}

std::vector<double> readSinusX(const char *fileName)
{
  std::ifstream ifs(fileName);
  if (!ifs)
  {
    std::ostringstream oss;
    oss << "readSinusX : Error while opening \"" << fileName << "\" file !";
    throw MyException(oss.str());
  }
  ifs.seekg(0, ifs.end);
  std::streampos length(ifs.tellg());
  ifs.seekg(0, ifs.beg);
  AutoPtr<char> data(new char[length]);
  std::vector<double> ret;
  while (!ifs.eof())
  {
    ifs.getline(data, length);
    std::streamsize szOfLine(ifs.gcount());
    double vals[3];
    if (isDataLine(data, szOfLine, vals[0], vals[1], vals[2]))
      ret.insert(ret.end(), vals, vals + 3);
  }
  return ret;
}

void performSubDiv(std::vector<double> &data, int nbOfSubDiv)
{
  constexpr int SPACEDIM = 3;
  if (nbOfSubDiv <= 1)
    return;
  if (data.size() % SPACEDIM != 0)
    throw MyException("Internal error : invalid size of data !");
  std::size_t nbPts(data.size() / SPACEDIM);
  if (nbPts <= 1)
    return;
  std::size_t newNbPts((nbOfSubDiv - 1) * (nbPts - 1) + nbPts);
  std::vector<double> newData(newNbPts * SPACEDIM);
  const double *inPt(data.data());
  double *pt(newData.data());
  for (auto i = 0; i < nbPts - 1; i++, inPt += SPACEDIM)
  {
    pt = std::copy(inPt, inPt + SPACEDIM, pt);
    const double *inPtNext(inPt + SPACEDIM);
    for (auto j = 1; j < nbOfSubDiv; j++)
    {
      double ratio((double)j / (double)nbOfSubDiv);
      pt = std::transform(inPt, inPt + SPACEDIM, inPtNext, pt, [ratio](const double &a, const double &b) { return a + ratio * (b - a); });
    }
  }
  pt = std::copy(inPt, inPt + SPACEDIM, pt);
  data = std::move(newData);
}

vtkStandardNewMacro(vtkSinusXReader);

vtkSinusXReader::vtkSinusXReader()
{
  this->SetNumberOfInputPorts(0);
}

int vtkSinusXReader::RequestInformation(vtkInformation *vtkNotUsed(request),
                                        vtkInformationVector **vtkNotUsed(inputVector),
                                        vtkInformationVector *outputVector)
{
  return 1;
}

int vtkSinusXReader::RequestData(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **vtkNotUsed(inputVector),
                                 vtkInformationVector *outputVector)
{
  vtkInformation *outInfo(outputVector->GetInformationObject(0));
  vtkPolyData *output(vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT())));
  //
  try
  {
    vtkNew<vtkPolyData> ret;
    if (!this->FileName)
      return 0;
    std::vector<double> data(readSinusX(this->FileName));
    performSubDiv(data, this->NumberOfSubdiv);
    if (data.size() % 3 != 0)
      throw MyException("Internal error : invalid size of data !");
    std::size_t nbPts(data.size() / 3);
    vtkNew<vtkDoubleArray> arr;
    arr->SetNumberOfComponents(3);
    arr->SetNumberOfTuples(nbPts);
    std::copy(data.begin(), data.end(), arr->GetPointer(0));
    vtkNew<vtkPoints> pts;
    pts->SetData(arr);
    ret->SetPoints(pts);
    vtkNew<vtkCellArray> verts;
    {
      vtkNew<vtkIdTypeArray> conn;
      conn->SetNumberOfComponents(1);
      conn->SetNumberOfTuples(2 * nbPts);
      vtkIdType *pt(conn->GetPointer(0));
      for (vtkIdType i = 0; i < nbPts; i++)
      {
        pt[2 * i] = 1;
        pt[2 * i + 1] = i;
      }
      verts->SetCells(nbPts, conn);
    }
    ret->SetVerts(verts);
    if (nbPts >= 1)
    {
      vtkNew<vtkCellArray> lines;
      {
        vtkNew<vtkIdTypeArray> conn;
        conn->SetNumberOfComponents(1);
        conn->SetNumberOfTuples(3 * (nbPts - 1));
        vtkIdType *pt(conn->GetPointer(0));
        for (vtkIdType i = 0; i < nbPts - 1; i++)
        {
          pt[3 * i] = 2;
          pt[3 * i + 1] = i;
          pt[3 * i + 2] = i + 1;
        }
        lines->SetCells(nbPts - 1, conn);
      }
      ret->SetLines(lines);
    }
    output->ShallowCopy(ret);
  }
  catch (MyException &e)
  {
    vtkErrorMacro(<< "vtkSinusXReader::RequestData : during read of " << this->FileName << " : " << e.what());
    return 0;
  }
  return 1;
}

void vtkSinusXReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
