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

/*=========================================================================

  Program:   ParaView
  Module:    vtkXYChartRepresentationColumns.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkXYChartRepresentationColumns
 *
 * vtkXYChartRepresentationColumns is a specialisation of vtkXYChartRepresentation
 * that supports column selection.
 */

#ifndef vtkXYChartRepresentationColumns_h
#define vtkXYChartRepresentationColumns_h

#include <vtkXYChartRepresentation.h>

class VTK_EXPORT vtkXYChartRepresentationColumns : public vtkXYChartRepresentation
{
public:
  static vtkXYChartRepresentationColumns* New();
  vtkTypeMacro(vtkXYChartRepresentationColumns, vtkXYChartRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  vtkXYChartRepresentationColumns() = default;
  ~vtkXYChartRepresentationColumns() override = default;

  void PrepareForRendering() override;

private:
  vtkXYChartRepresentationColumns(const vtkXYChartRepresentationColumns&) = delete;
  void operator=(const vtkXYChartRepresentationColumns&) = delete;
};

#endif // vtkXYChartRepresentationColumns_h
