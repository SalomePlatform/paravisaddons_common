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

#ifndef __vtkPVMetaDataInformationERIO_h
#define __vtkPVMetaDataInformationERIO_h

#include "vtkPVInformation.h"

class vtkDataObject;
class vtkInformationDataObjectKey;

class VTK_EXPORT vtkPVMetaDataInformationERIO : public vtkPVInformation
{
public:
  static vtkPVMetaDataInformationERIO* New();
  vtkTypeMacro(vtkPVMetaDataInformationERIO, vtkPVInformation);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Transfer information about a single object into this object.
  virtual void CopyFromObject(vtkObject*);

  //BTX
  // Description:
  // Manage a serialized version of the information.
  virtual void CopyToStream(vtkClientServerStream*);
  virtual void CopyFromStream(const vtkClientServerStream*);
  virtual void AddInformation(vtkPVInformation*);
  //ETX

  // Description:
  // Returns the Information Data.
  vtkGetObjectMacro(InformationData, vtkDataObject);

//BTX
protected:
  vtkPVMetaDataInformationERIO();
  ~vtkPVMetaDataInformationERIO();
  void SetInformationData(vtkDataObject*);
  vtkDataObject* InformationData;

private:
  vtkPVMetaDataInformationERIO(const vtkPVMetaDataInformationERIO&); // Not implemented
  void operator=(const vtkPVMetaDataInformationERIO&); // Not implemented
//ETX
};

#endif
