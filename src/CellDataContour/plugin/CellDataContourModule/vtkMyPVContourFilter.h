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

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMyContourFilter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkMyPVContourFilter
 * @brief   generate isosurfaces/isolines from scalar values
 *
 * vtkMyPVContourFilter is an extension to vtkMyContourFilter. It adds the
 * ability to generate isosurfaces / isolines for AMR dataset.
 *
 * @warning
 * Certain flags in vtkAMRDualContour are assumed to be ON.
 *
 * @sa
 * vtkMyContourFilter vtkAMRDualContour
*/

#ifndef vtkMyPVContourFilter_h
#define vtkMyPVContourFilter_h

#include "vtkMyContourFilter.h"

class VTK_EXPORT vtkMyPVContourFilter : public vtkMyContourFilter
{
public:
  vtkTypeMacro(vtkMyPVContourFilter, vtkMyContourFilter);

  void PrintSelf(ostream &os, vtkIndent indent) override;

  static vtkMyPVContourFilter* New();

  int ProcessRequest(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

protected:
  vtkMyPVContourFilter();
  ~vtkMyPVContourFilter() override;

  int RequestData(vtkInformation *request, vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  virtual int RequestDataObject(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector);

  int FillInputPortInformation(int port, vtkInformation* info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * Class superclass request data. Also handles iterating over
   * vtkHierarchicalBoxDataSet.
   */
  int ContourUsingSuperclass(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector);

private:
  vtkMyPVContourFilter(const vtkMyPVContourFilter&) = delete;
  void operator=(const vtkMyPVContourFilter&) = delete;
};

#endif // vtkMyPVContourFilter_h
