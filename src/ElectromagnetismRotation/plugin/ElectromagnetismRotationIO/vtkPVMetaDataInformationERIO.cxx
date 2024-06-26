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
// Author : Anthony Geay

#include "vtkPVMetaDataInformationERIO.h"

#include "vtkAlgorithm.h"
#include "vtkAlgorithmOutput.h"
#include "vtkClientServerStream.h"
#include "vtkExecutive.h"
#include "vtkDataObject.h"
#include "vtkGenericDataObjectReader.h"
#include "vtkGenericDataObjectWriter.h"
#include "vtkInformationDataObjectMetaDataKey.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"

#include "MEDCouplingRefCountObject.hxx"

#include <sstream>

vtkStandardNewMacro(vtkPVMetaDataInformationERIO);
vtkCxxSetObjectMacro(vtkPVMetaDataInformationERIO, InformationData, vtkDataObject);

static vtkInformationDataObjectMetaDataKey* GetMEDReaderMetaDataIfAny()
{
  static const char ZE_KEY[] = "vtkFileSeriesGroupReader::META_DATA";
  MEDCoupling::GlobalDict* gd(MEDCoupling::GlobalDict::GetInstance());
  if (!gd->hasKey(ZE_KEY))
    return 0;
  std::string ptSt(gd->value(ZE_KEY));
  void* pt(0);
  std::istringstream iss(ptSt);
  iss >> pt;
  return reinterpret_cast<vtkInformationDataObjectMetaDataKey*>(pt);
}

//----------------------------------------------------------------------------
vtkPVMetaDataInformationERIO::vtkPVMetaDataInformationERIO()
{
  this->InformationData = NULL;
}

//----------------------------------------------------------------------------
vtkPVMetaDataInformationERIO::~vtkPVMetaDataInformationERIO()
{
  this->SetInformationData(NULL);
}

//----------------------------------------------------------------------------
void vtkPVMetaDataInformationERIO::CopyFromObject(vtkObject* obj)
{
  this->SetInformationData(NULL);

  vtkAlgorithmOutput* algOutput = vtkAlgorithmOutput::SafeDownCast(obj);
  if (!algOutput)
    {
    vtkAlgorithm* alg = vtkAlgorithm::SafeDownCast(obj);
    if (alg)
      {
      algOutput = alg->GetOutputPort(0);
      }

    }
  if (!algOutput)
    {
    vtkErrorMacro("Information can only be gathered from a vtkAlgorithmOutput.");
    return;
    }

  vtkAlgorithm* reader = algOutput->GetProducer();
  vtkInformation* info = reader->GetExecutive()->GetOutputInformation(
    algOutput->GetIndex());

  if (info && info->Has(GetMEDReaderMetaDataIfAny()))
    {
    this->SetInformationData(vtkDataObject::SafeDownCast(info->Get(GetMEDReaderMetaDataIfAny())));
    }
}

//----------------------------------------------------------------------------
void vtkPVMetaDataInformationERIO::CopyToStream(vtkClientServerStream* css)
{
  css->Reset();
  if (!this->InformationData)
    {
    *css << vtkClientServerStream::Reply
         << vtkClientServerStream::InsertArray(
           static_cast<unsigned char*>(NULL), 0)
         << vtkClientServerStream::End;
    return;
    }

  vtkDataObject* clone = this->InformationData->NewInstance();
  clone->ShallowCopy(this->InformationData);

  vtkGenericDataObjectWriter* writer = vtkGenericDataObjectWriter::New();
  writer->SetFileTypeToBinary();
  writer->WriteToOutputStringOn();
  writer->SetInputData(clone);
  writer->Write();

  *css << vtkClientServerStream::Reply
       << vtkClientServerStream::InsertArray(
         writer->GetBinaryOutputString(),
         writer->GetOutputStringLength())
       << vtkClientServerStream::End;
  writer->RemoveAllInputs();
  writer->Delete();
  clone->Delete();
}

//----------------------------------------------------------------------------
void vtkPVMetaDataInformationERIO::CopyFromStream(const vtkClientServerStream* css)
{
  this->SetInformationData(0);
  vtkTypeUInt32 length;
  if (css->GetArgumentLength(0, 0, &length) && length > 0)
    {
    unsigned char* raw_data = new unsigned char[length];
    css->GetArgument(0, 0, raw_data, length);
    vtkGenericDataObjectReader* reader = vtkGenericDataObjectReader::New();
    reader->SetBinaryInputString(reinterpret_cast<const char*>(raw_data), length);
    reader->ReadFromInputStringOn();
    delete []raw_data;
    reader->Update();
    this->SetInformationData(reader->GetOutput());
    reader->Delete();
    }
}

void vtkPVMetaDataInformationERIO::AddInformation(vtkPVInformation*)
{
}

//----------------------------------------------------------------------------
void vtkPVMetaDataInformationERIO::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "InformationData: " <<  this->InformationData << endl;
}
