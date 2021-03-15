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

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSerafinReader.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

////// Reader for files 2D/3D generated by The  TELEMAC modelling system \\\\\
// Module developped by herve ozdoba - Sept 2008 ( herve-externe.ozdoba at edf.fr / herve at ozdoba.fr )
// Please address all comments to Regina Nebauer ( regina.nebauer at edf.fr )
// >>> Test version

#include "FFileReader.h"
#include "stdSerafinReader.h"
#include "vtkSerafinReader.h"

#include "vtkErrorCode.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkCellArray.h"
#include "vtkNew.h"

#include <iostream>
#include <list>
#ifdef WIN32
#define NOMINMAX
#endif

class vtkSerafinReader::vtkInternal
{
public:
  vtkInternal() { }
  vtkUnstructuredGrid *getPointer() const { return _usg.GetPointer(); }
  void geometryHasBeenRead() { _geometry_read=true; }
  bool hasGeometryAlreadyRead() const { return _geometry_read; }
private:
  bool _geometry_read = false;
  vtkNew<vtkUnstructuredGrid> _usg;
};

/** +++++++++++++++++ Définition des méthodes de la classe FFileReader +++++++++++++++++ **/

/* ******************* Constructeur *******************
 * Ce constructeur reçoit un flux de lecture de fichier en argument .
 * Gloabalement, les initialisations sont effectuée ici et le premier entier est lu pour
 * déterminer dans quelle configuration d'écriture on se place en association avec le fichier .
 * La différenciation petit/grand boutien est faite à la lecture du premier entier
 * qui doit valoir la taille maximale du titre soit 80 caractères (à l'heure où j'écris ces lignes) .
 */
FFileReader :: FFileReader(ifstream* stream)
{
  // différentes initialisations
  this->BigEndian  = false ;
  this->BlocSize   = 0 ;
  this->FileStream = stream;

  // lecture de l'entête
  readBlocSize ();
  //vtkDebugMacro( << "BlocSize Little Endian " << BlocSize << "\n");
  if (this->BlocSize != TITLE_MAX_SIZE) // pas d'échange d'octet à la lecture
  {
    this->BigEndian = true ;
    // Relecture de l'entête
    readBlocSize ();
    //vtkDebugMacro( << "BlocSize Big Endian " << BlocSize << "\n");
  }

  // The float reader will be set later when we identify if single or double
  FFileReader::readIntArray   = &FFileReader::g_readInt32Array   ;
};

/* Lecture d'un tableau d'entier arr de taille size avec inversion de octets */
int64_t FFileReader :: s_readInt32Array(int* arr, const int size)
{
  FileStream->read ((char*)(arr), sizeof(int)*size);
  Swap32Array (size, (char*)(arr));
  return FileStream->tellg();
};

/* Lecture d'un tableau d'entier arr de taille size sans inversion de octets */
int64_t FFileReader :: ns_readInt32Array(int* arr, const int size)
{
  FileStream->read ((char*)(arr), sizeof(int)*size);
  return FileStream->tellg();
};

/* Lecture d'un tableau de flottants arr de taille size avec inversion de octets */
int64_t FFileReader :: s_readFloat32Array(double* arr, const int size)
{
  float *tmp = new float[size];
  FileStream->read ((char*)(tmp), sizeof(float)*size);
  Swap32Array (size, (char*)(tmp));
  for(int i; i<size; i++) arr[i] = (double) tmp[i];
  return FileStream->tellg();
};

/* Lecture d'un tableau de flottants arr de taille size sans inversion de octets */
int64_t FFileReader :: ns_readFloat32Array(double* arr, const int size)
{
  float *tmp = new float[size];
  FileStream->read ((char*)(tmp), sizeof(float)*size);
  for(int i; i<size; i++) arr[i] = (double) tmp[i];
  return FileStream->tellg();
};
// generic single entry read functions
int64_t FFileReader :: g_readInt32Array(int* arr, const int size)
{
  if ( this->BigEndian ) {
    return s_readInt32Array(arr, size);
  }
  return ns_readInt32Array(arr, size);
};

int64_t FFileReader :: g_readFloat32Array(double* arr, const int size)
{
  if ( this->BigEndian ) {
    return s_readFloat32Array(arr, size);
  }
  return ns_readFloat32Array(arr, size);
};

/* Lecture d'un tableau d'entier arr de taille size avec inversion de octets */
int64_t FFileReader :: s_readInt64Array(int64_t* arr, const int size)
{
  FileStream->read ((char*)(arr), sizeof(int64_t)*size);
  Swap64Array (size, (char*)(arr));
  return FileStream->tellg();
};

/* Lecture d'un tableau d'entier arr de taille size sans inversion de octets */
int64_t FFileReader :: ns_readInt64Array(int64_t* arr, const int size)
{
  FileStream->read ((char*)(arr), sizeof(int64_t)*size);
  return FileStream->tellg();
};

/* Lecture d'un tableau de flottants arr de taille size avec inversion de octets */
int64_t FFileReader :: s_readFloat64Array(double* arr, const int size)
{
  FileStream->read ((char*)(arr), sizeof(double)*size);
  Swap64Array (size, (char*)(arr));
  return FileStream->tellg();
};

/* Lecture d'un tableau de flottants arr de taille size sans inversion de octets */
int64_t FFileReader :: ns_readFloat64Array(double* arr, const int size)
{
  FileStream->read ((char*)(arr), sizeof(double)*size);
  return FileStream->tellg();
};
// generic single entry read functions
int64_t FFileReader :: g_readInt64Array(int64_t* arr, const int size)
{
  if ( this->BigEndian ) {
    return s_readInt64Array(arr, size);
  }
  return ns_readInt64Array(arr, size);
};

int64_t FFileReader :: g_readFloat64Array(double* arr, const int size)
{
  if ( this->BigEndian ) {
    return s_readFloat64Array(arr, size);
  }
  return ns_readFloat64Array(arr, size);
};

/* ******************* Destructeur ***************** */
// TODO compléter cette méthode !!!
FFileReader :: ~FFileReader()
{
  // Ne rien faire pour le moment
};

/** +++++++++++++++++ Définition des méthodes de la classe stdSerafinReader +++++++++++++++++ **/

/* ******************* Constructeur ***************** */
stdSerafinReader :: stdSerafinReader(ifstream* stream, int BuildVectors) : FFileReader(stream)
{
  // TODO Initialisation des variables
  this->metadata  = new SerafinMetaData();
  this->index  = new SerafinIndexInfo();

  // Lecture des metadonnée
  //vtkDebugMacro( << "Reafin metadata" << endl);
  this->readMetaData ();

  //Création de l'index
  //vtkDebugMacro( << "Creating Index" << endl);
  this->createIndex ();

  // Identify variable vector
  //vtkDebugMacro( << "Computing Var Infor" << endl);
  this->ComputeVarInfo(BuildVectors);
};

/* ******************* Destructeur ***************** */
// TODO compléter cette méthode !!!
//stdSerafinReader :: ~stdSerafinReader()
//{
//  // Ne rien faire pour le moment, provoque une 'legere fuite memoire' maitrisee
//};
void stdSerafinReader::ComputeVarInfo(int BuildVectors)
{
  int pos, ncomp;
  int found;

  if (Is3Dfile())
    ncomp = 3;
  else
    ncomp = 2;

  for( int i; i < this->metadata->VarNumber ; i++) {
    // Only converting to vector if asked for
    if(! BuildVectors){
      metadata->nVarList[i].ncomp = 0;
      metadata->nVarList[i].icomp = 0;
      continue;
    }
    string name(metadata->nVarList[i].name);

    // Exception for which we do not create a vector
    list<string> exceptions {"COTE Z", "ELEVATION Z"};

    found = 0;
    for(auto const &excep : exceptions){
      pos = name.find(excep);
      if (pos != std::string::npos){
        metadata->nVarList[i].ncomp = 0;
        metadata->nVarList[i].icomp = 0;
        found = 1;
        continue;
      }
    }
    if (found) continue;

    // List of ends of variable that should be converted to vectors
    list<string> vec0 {" U ", " X ", "QX ", "U0 "};
    list<string> vec1 {" V ", " Y ", "QY ", "V0 "};
    list<string> vec2 {" W ", " Z ", "QZ ", "W0 "};
    list<list<string>> vec {vec0, vec1, vec2};

    found = 0;
    int k = 0;
    for(auto const &veci : vec){
      for(auto const &str : veci){
        pos = name.find(str);
        if (pos != std::string::npos){
          metadata->nVarList[i].ncomp = ncomp;
          metadata->nVarList[i].icomp = k;
          metadata->nVarList[i].name[pos+1] = '*';
          if (str[0] != ' ') metadata->nVarList[i].name[pos] = '*';
          found = 1;
          break;
        };
      }
      k++;
      if (found) break;
    }
    // Found a vector go to next variable
    if (found) continue;
    // Default values
    metadata->nVarList[i].ncomp = 0;
    metadata->nVarList[i].icomp = 0;
  }
};

/* ******************* createIndex ***************** */
/* Cette méthode crée un index de taille et de position à partir des informations meta
 * afin de faciliter la lecture du fichier serafin .
 */
void stdSerafinReader::createIndex ()
{
  int tag = 0 ;
  int nnodes = GetNumberOfNodes();
  int ndp = GetNodeByElements();
  int nelem = GetNumberOfElement();

  // TODO: Identify FloatSize (read tag of coordinates)
  this->index->IntSize = sizeof(int) ;
  this->index->TagSize = sizeof(int) ;

  this->index->FileSize = GetFileSize();
  this->index->MetaSize = FileStream->tellg();

  this->index->ConnectivityPosition  = this->index->MetaSize;

  this->index->XPosition = (this->index->MetaSize) +
                           (this->index->IntSize*nnodes+2*this->index->TagSize) +
                           (this->index->IntSize*ndp*nelem+2*this->index->TagSize);

  // Identifying float precision from tag of coordinates
  FileStream->seekg( this->index->XPosition);
  (*this.*readIntArray)(&tag, 1);
  this->index->FloatSize = int(tag/nnodes) ;
  //vtkDebugMacro( << "Float Size: " << this->index->FloatSize << endl);
  if (this->index->FloatSize == 4){
    FFileReader::readFloatArray = &FFileReader::g_readFloat32Array ;
  }else{
    FFileReader::readFloatArray = &FFileReader::g_readFloat64Array ;
  }

  // Size of a time info
  this->index->TimeSize     = 2*this->index->TagSize + this->index->FloatSize;
  // Size of a Field
  this->index->FieldSize    = this->index->FloatSize*nnodes+2*this->index->TagSize;

  // Size of the whole data
  this->index->DataSize     = (this->index->FileSize)
                              - (index->MetaSize)
                              - (this->index->IntSize*nnodes+2*this->index->TagSize)
                              - 2*index->FieldSize
                              - (this->index->IntSize*ndp*nelem+2*this->index->TagSize);
  // Index to data
  this->index->DataPosision = (this->index->FileSize) - (this->index->DataSize);

  // Index of Y coordinates
  this->index->YPosition = this->index->XPosition + this->index->FieldSize;
  // Size of data bloc for one time step
  this->index->DataBlocSize = this->index->TimeSize + GetNumberOfVars()*this->index->FieldSize;

  /*............................................................................................*/


  this->index->NumberOfDate     = (this->index->DataSize)/(this->index->DataBlocSize);
  //vtkDebugMacro(<< "Number of Date: " << this->index->NumberOfDate << endl);

};

/* ******************* readMetaData ***************** */
/* Cette methode permet de de lire les métadata dans le but de recueillir les informations
 * essentielles incluses dans le fichier . Globalement, la demarche sequentielle est la suivante :
 *   -  lecture du titre et suppression des espaces en fin de chaine s'il y en a
 *   -  lecture du nombre de variables
 *   -  lecture du nom des variables et des leurs unités respectives
 *   -  lecture des paramamètres
 *  -  lecture des informations de discrétisation
 */
int stdSerafinReader :: readMetaData ()
{

  //Lecture du titre
  if (ReadString(metadata->Title, TITLE_MAX_SIZE) != 88) return 0;// metadata->Title[TITLE_MAX_SIZE]='\0';
  DeleteBlank(metadata->Title, TITLE_MAX_SIZE-8);

  //lecture du nombre de variables (on passe les entete)
  skipReadingHeader(FileStream);   //skip reclen
  // read linear varsno
  if ((*this.*readIntArray)(&(metadata->VarNumber), 1) != 96) return 0;
  skipReadingHeader(FileStream);  // skip quad varno
  skipReadingHeader(FileStream);  // skip reclen

  //lecture des variables
  metadata->VarList = (char *)new SerafinVar[metadata->VarNumber];
  metadata->nVarList = new SerafinVar[metadata->VarNumber];

  //vtkDebugMacro( << "nVarList Size " << sizeof(metadata->nVarList) << "\n");

  {
    int compteur = 0 ;
    char buffer[VAR_DESC_SIZE*2];
    for( compteur; compteur < metadata->VarNumber ; compteur++) {
      // must read full buffer as each varaible is a file record
      ReadString(&buffer[0], VAR_DESC_SIZE*2);
      strncpy(metadata->nVarList[compteur].name, &buffer[0], VAR_DESC_SIZE);
      strncpy(metadata->nVarList[compteur].unit, &buffer[VAR_DESC_SIZE], VAR_DESC_SIZE);
      metadata->nVarList[compteur].name[VAR_DESC_SIZE]='\0';
      metadata->nVarList[compteur].unit[VAR_DESC_SIZE]='\0';
    }
  };

  // Lecture des parametres et, si necessaire, de la date de simu
  skipReadingHeader(FileStream);
  (*this.*readIntArray)(metadata->IParam, PARAM_NUMBER);
  skipReadingHeader(FileStream);

  if (metadata->IParam[9] == 1)// Si la date est indiquée
  {
    skipReadingHeader(FileStream);
    (*this.*readIntArray)(metadata->Date, DATE_NUMBER);
    skipReadingHeader(FileStream);
  };

  //lecture des information de discrietisation
  skipReadingHeader(FileStream);
  (*this.*readIntArray)(metadata->DiscretizationInfo, DISC_DESC_SIZE);
  skipReadingHeader(FileStream);

  // On lit l'entete du bloc de lecture pour connaitre la taille de la table de connectivite
  if (IsBigEndian()) s_readBlocSize ();else  ns_readBlocSize ();


  return FileStream->tellg();
};

/** +++++++++++++++++ Définition des méthodes de la classe vtkSerafinReader +++++++++++++++++ **/

#include "vtkObjectFactory.h"

//vtkCxxRevisionMacro(vtkSerafinReader, "$Revision: 0.2 $");
vtkStandardNewMacro(vtkSerafinReader);

vtkSerafinReader::vtkSerafinReader():Internal(nullptr)
{

  //vtkDebugMacro( << "Instanciation du lecteur Serafin");

  this->FileName    = NULL;
  this->FileStream  = NULL;
  this->Reader      = NULL;
  this->TimeStep    = 0;
  this->Internal=new vtkInternal;

  this->SetNumberOfInputPorts(0);
};

vtkSerafinReader::~vtkSerafinReader()
{
  if (this->FileName)
  {
    this->SetFileName(0);
  }
  delete this->Internal;

}

int vtkSerafinReader::RequestInformation(vtkInformation *vtkNotUsed(request),
           vtkInformationVector **vtkNotUsed(inputVector),
           vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  //outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),1);

  if ( !this->FileName )
  {
    vtkErrorMacro("No filename specified");
    return 0;
  }

  this->FileStream = new ifstream(this->FileName, ifstream::binary|ifstream::in);

  if (this->FileStream->fail())
  {
    this->SetErrorCode(vtkErrorCode::FileNotFoundError);
    delete this->FileStream;
    this->FileStream = NULL;
    vtkErrorMacro("Specified filename not found");
    return 0;
  }

  this->Reader   = new stdSerafinReader( FileStream, this->BuildVectors);

  {//Gestion du temps
    const int totime = this->Reader->GetTotalTime();
    if (totime > 1)
    {
      int i=0;
      double* TimeValues = new double[totime];

      for (i=0; i<totime ;i++) {TimeValues[i] = this->Reader->GetTime(i) ;}

      outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &TimeValues[0],  totime);

      double timeRange[2];
      timeRange[0] = TimeValues[0];
      timeRange[1] = TimeValues[totime-1];
      outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),  timeRange, 2);

    };
  }
  return 1;
}

int vtkSerafinReader::RequestData(vtkInformation *vtkNotUsed(request),
          vtkInformationVector **vtkNotUsed(inputVector),
          vtkInformationVector *outputVector)
{
  int totime = this->Reader->GetTotalTime();
  vtkInformation     *outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid   *output = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  int tsLength = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  double *steps =   outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  double requestedTimeSteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

  if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) && tsLength>0)
  {
    // Get the requested time step. We only support requests of a single time
    // step in this reader right now
    double requestedTimeSteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

    // find the first time value larger than requested time value
    // this logic could be improved
    int cnt = 0;
    while (cnt < tsLength-1 && steps[cnt] < requestedTimeSteps)
    {
      cnt++;
    }

    this->TimeStep = cnt;
  }

  //vtkDebugMacro( << "Serafin steps <" << steps << ">..." << requestedTimeSteps << this->TimeStep);

  if ( outInfo->Has( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() ) )
  {
    double* steps = outInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
  };


  if ( this->FileStream == NULL )
  {
    return 0;
  }

  //vtkDebugMacro( << "Reading Time " << this->TimeStep << endl);

  // Lecture de la geometrie
  if(!this->Internal->hasGeometryAlreadyRead())
    this->ReadGeometry(Internal->getPointer(), this->TimeStep);
  Internal->geometryHasBeenRead();
  output->ShallowCopy(Internal->getPointer());

  // Lecture des donnees
  //vtkDebugMacro( << "Reading Data " << this->TimeStep << endl);
  this->ReadData(output,  this->TimeStep);
  //vtkDebugMacro( << "Request done" << endl);

  return 1;
}

void vtkSerafinReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: "           << (this->FileName ? this->FileName : "(none)") << endl;
  os << indent << "Number Of Nodes: "     << this->Reader->GetNumberOfNodes()    << endl;
  os << indent << "Number Of Node Fields: "     << this->Reader->GetNumberOfVars()     << endl;
  os << indent << "Number Of Cells: "     << this->Reader->GetNumberOfElement()     << endl;
}

void vtkSerafinReader::ReadGeometry(vtkUnstructuredGrid *output, int time)
{
  vtkDoubleArray *coords = vtkDoubleArray::New();
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(this->Reader->GetNumberOfNodes());

  //vtkDebugMacro( << "Reading coordinates" << endl);
  this->Reader->WriteCoord(coords->GetPointer (0), time);

  //Lecture de la table de connectivite
  //vtkDebugMacro( << "Reading connectivity" << endl);
  {
    int i = 0, k = 0, l = 0;
    vtkIdType list[27];
    const int size = this->Reader->GetNodeByElements()*this->Reader->GetNumberOfElement();
    int* arr = new int[size];

    switch(this->Reader->GetNodeByElements())
    {
      case 3 : l = VTK_TRIANGLE ; break;
      case 4 : l = (this->Reader->Is3Dfile())? VTK_TETRA : VTK_QUAD ;  break;
      case 5 : l = VTK_PYRAMID; break;
      case 6 : l = VTK_WEDGE;  break;
      case 8 : l = VTK_HEXAHEDRON ;break;
      default:
      {
        vtkErrorMacro( << "cell type is not supported\n");return;
      }
    }

    //vtkDebugMacro( << "Writting connectivity\n");
    this->Reader->WriteConnectivity(arr);
    output->Allocate(this->Reader->GetNumberOfNodes(), this->Reader->GetNumberOfNodes());

    for(i = 0; i < this->Reader->GetNumberOfElement(); i++)
    {
      for(k = 0; k < this->Reader->GetNodeByElements(); k++)
        list[k] = arr[this->Reader->GetNodeByElements()*i+k]-1;

      output->InsertNextCell(l, this->Reader->GetNodeByElements(), list);
    };

    delete[] arr;
  };

  //vtkDebugMacro( << "Setting points\n");
  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  coords->Delete();

  output->SetPoints(points);
  points->Delete();
  //vtkDebugMacro( << "Read Geometry done\n");

}



void vtkSerafinReader::ReadData(vtkUnstructuredGrid *output, int time)
{
  int i = 0, dim = 1;int vel =0 ;
  char name[VAR_DESC_SIZE+1];

  const int size = this->Reader->GetNumberOfNodes();

  const int ideb = 0;
  const int ifin = this->Reader->GetNumberOfVars();
  SerafinVar * var ;

  for (i = ideb ; i<ifin ; i++)
  {
    var = &(this->Reader->metadata->nVarList[i]);
    //vtkDebugMacro( << "ReadData varname" << endl);
    //vtkDebugMacro( << "  id: " << i << endl);
    //vtkDebugMacro( << "  name: *" << var->name << "*\n");
    //vtkDebugMacro( << "  icomp/ncomp: " << var->icomp << "/" << var->ncomp << endl);
    vtkDoubleArray *data = vtkDoubleArray::New();
    std::string name(var->name);
    name = name.substr(0, name.find_last_not_of(" \n")+1);

    if (var->ncomp != 0){
      data->SetName(name.c_str());
      data->SetNumberOfComponents(3);
      data->SetNumberOfTuples(size);
      this->Reader->GetVarRangeValues(size, var->ncomp, i, data->GetPointer(0), time);
      i+= var->ncomp-1;

    }else{
      data->SetName(name.c_str());
      data->SetNumberOfComponents(1);
      data->SetNumberOfTuples(size);
      this->Reader->GetVarValues(time, i, 0, data->GetPointer (0), size);
    }

    output->GetPointData()->AddArray(data);

    data->Delete();

  }
}
