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

/*=========================================================================

  Program:   ParaView
  Module:    vtkAppendAttributesOverTime.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkAppendAttributesOverTime
 * @brief   Multiple inputs with one output.
 *
 * vtkAppendAttributesOverTime put all arrays from all inputs into one output,
 * and this for each timesteps.
 * The output data object is the same as the first data input, and has no timeteps.
 * It extends vtkMergeArrays to process every available timestep.
 * The new arrays will have the name mangled to be the original array name plus
 * `_input_<inputid>_ts_<tsindex>` where `<inputid>` is the id/index of the input filter
 * that is providing that array, and `<tsindex>` is the timestep index of corresponing input.
 */

#ifndef vtkAppendAttributesOverTime_h
#define vtkAppendAttributesOverTime_h

#include <vtkMergeArrays.h>

#include <vtkInformation.h>
#include <vtkNew.h>                          // vtkNew
#include <vtkSmartPointer.h>                 // Smart Pointer

#include <vector>

class vtkFieldData;
class vtkDataObject;

class VTK_EXPORT vtkAppendAttributesOverTime : public vtkMergeArrays
{
public:
  static vtkAppendAttributesOverTime* New();
  vtkTypeMacro(vtkAppendAttributesOverTime, vtkMergeArrays);

protected:
  vtkAppendAttributesOverTime();
  ~vtkAppendAttributesOverTime() override = default;

  //@{
  /**
   * Given an array name, return an appropriate name to use for the output array.
   * Reimplemented to add original timestep index in the name.
   * Returns true.
   */
  bool GetOutputArrayName(vtkFieldData* arrays, const char* inArrayName, int inputIndex,
    std::string& outArrayName) override;
  //@}

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  vtkAppendAttributesOverTime(const vtkAppendAttributesOverTime&) = delete;
  void operator=(const vtkAppendAttributesOverTime&) = delete;

  int UpdateTimeIndex = 0;
  int CurrentInputIndex = 0;

  /**
   * this->TimeSteps contains for each input a vector of timesteps
   * time = this->TimeSteps[inputIndex][tsIndex]
   */
  std::vector<std::vector<double> > TimeSteps;
  vtkNew<vtkInformation> CurrentOutInfo;
  vtkSmartPointer<vtkDataObject> TempDataObject;
  double OriginalTimeStep;
  bool RestoreOriginalTimeStep;
};

#endif // vtkAppendAttributesOverTime_h
