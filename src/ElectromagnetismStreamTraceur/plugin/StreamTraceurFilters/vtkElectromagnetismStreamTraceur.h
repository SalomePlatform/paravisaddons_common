// Copyright (C) 2021-2026  CEA, EDF
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

#pragma once

#include "vtkFiltersFlowPathsModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

#include "vtkInitialValueProblemSolver.h" // Needed for constants

class vtkAbstractInterpolatedVelocityField;
class vtkCompositeDataSet;
class vtkDataArray;
class vtkDataSetAttributes;
class vtkDoubleArray;
class vtkExecutive;
class vtkGenericCell;
class vtkIdList;
class vtkIntArray;
class vtkPoints;

#include <vector>

class VTK_EXPORT vtkElectromagnetismStreamTraceur : public vtkPolyDataAlgorithm
{
public:
  enum
  {
    FORWARD,
    BACKWARD,
    BOTH
  };

    enum Solvers
  {
    RUNGE_KUTTA2,
    RUNGE_KUTTA4,
    RUNGE_KUTTA45,
    NONE,
    UNKNOWN
  };


  vtkTypeMacro(vtkElectromagnetismStreamTraceur, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  static vtkElectromagnetismStreamTraceur* New();


  void SetSourceData(vtkDataSet* source);
  vtkDataSet* GetSource();
  void SetSourceConnection(vtkAlgorithmOutput* algOutput);

  vtkSetClampMacro(IntegrationDirection, int, FORWARD, BOTH);
  vtkGetMacro(IntegrationDirection, int);

  void SetIntegratorType(int type) { IntegratorType=type; }
  
  void SetIntegrationStepUnit(int unit) { IntegrationStepUnit=unit; }
  
  vtkSetMacro(InitialIntegrationStep, double);
  vtkGetMacro(InitialIntegrationStep, double);

  vtkSetMacro(MinimumIntegrationStep, double);
  vtkGetMacro(MinimumIntegrationStep, double);

  vtkSetMacro(MaximumIntegrationStep, double);
  vtkGetMacro(MaximumIntegrationStep, double);

  vtkSetMacro(TerminalSpeed, double);
  vtkGetMacro(TerminalSpeed, double);

  vtkSetMacro(MaximumError, double);
  vtkGetMacro(MaximumError, double);

  vtkSetMacro(MaximumNumberOfSteps, vtkIdType);
  vtkGetMacro(MaximumNumberOfSteps, vtkIdType);

protected:
  vtkElectromagnetismStreamTraceur();
  ~vtkElectromagnetismStreamTraceur() override = default;
  // Create a default executive.
  vtkExecutive* CreateDefaultExecutive() override;
  static const char *GetColorArrayName();
  // hide the superclass' AddInput() from the user and the compiler
  void AddInput(vtkDataObject*) { vtkErrorMacro(<< "AddInput() must be called with a vtkDataSet not a vtkDataObject."); }

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int FillInputPortInformation(int, vtkInformation*) override;

  int SetupOutput(vtkInformation* inInfo, vtkInformation* outInfo);

  vtkCompositeDataSet* InputData;
private:
  vtkElectromagnetismStreamTraceur(const vtkElectromagnetismStreamTraceur&) = delete;
  void operator=(const vtkElectromagnetismStreamTraceur&) = delete;
  static const char NAME_COLOR_ARRAY[];
public:
  int IntegrationDirection;
  int IntegratorType;
  int IntegrationStepUnit;
  double InitialIntegrationStep; 
  double MinimumIntegrationStep;
  double MaximumIntegrationStep;
  double TerminalSpeed;
  double MaximumError;
  vtkIdType MaximumNumberOfSteps;
};
