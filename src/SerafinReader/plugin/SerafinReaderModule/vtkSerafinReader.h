// Copyright (C) 2021-2022  CEA/DEN, EDF R&D
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
  Module:    $RCSfile: vtkSerafinReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

////// Reader for files generated by The  TELEMAC modelling system \\\\\
// Module developped by herve ozdoba - Sept 2008 ( herve-externe.ozdoba at edf.fr / herve at ozdoba.fr )
// Please address all comments to Regina Nebauer ( regina.nebauer at edf.fr )
// >>> Test version

#ifndef __vtkSerafinReader_h__
#define __vtkSerafinReader_h__

/** -- Inclusions issues de la bibliotheque standard du C++ -- */

#include <fstream>
#include <string>
#include <cstdio>
#include <iostream>
#include <cstring>

using namespace std;

/** -- Inclusion des entetes de la bibliotheque vtk -- **/

#include "vtkUnstructuredGridAlgorithm.h"

#include "vtkStringArray.h"


#include "stdSerafinReader.h"

#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkStdString.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkCellArray.h"


/** ********************************************************************************************* **/
/** -- Definition de la classe de lecture des fichiers externes au format Serafin pour Telemac -- **/
/** ********************************************************************************************* **/

class VTK_EXPORT vtkSerafinReader  : public vtkUnstructuredGridAlgorithm
{
public:

  static vtkSerafinReader *New();

  vtkTypeMacro(vtkSerafinReader,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Mandatory for time manipulation
  vtkSetMacro(TimeStep, int);
  vtkGetMacro(TimeStep, int);
  int GetNumberOfTimeSteps() {return this->Reader->GetTotalTime();};

  vtkSetMacro(BuildVectors, int);
  vtkGetMacro(BuildVectors, int);

protected:

  // Implementation du constructeur associe a la classe
  vtkSerafinReader();

  // Implementation du descructeur
  ~vtkSerafinReader();

  int RequestInformation  (vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestData    (vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  // Lecture de la geometrie du maillage
  void ReadGeometry   (vtkUnstructuredGrid *output, int time);

  // Lecture des donnees de la simulation au niveau des noeuds et des cellules.
  void ReadData   (vtkUnstructuredGrid *output, int time);

  char    *FileName;   // Nom du fichier ouvert par le logiciel Paraview
  ifstream *FileStream;// Flux de lecture du fichier

  int TimeStep; // Internal variable to remember actual time step
  int BuildVectors; // Defines if variables are to be converter to vectors (when necessary)

  stdSerafinReader* Reader; /** /!\ Instance de lecture du fichier Serafin **/

  class vtkInternal;
  vtkInternal *Internal;

private:
  vtkSerafinReader(const vtkSerafinReader&);  // Pas implemente
  void operator=(const vtkSerafinReader&);    // Pas implemente

}; /* class_vtkSerafinReader */

#endif
