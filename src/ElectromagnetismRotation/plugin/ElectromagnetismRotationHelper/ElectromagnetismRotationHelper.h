// Copyright (C) 2021-2023  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#pragma once

#include <sstream>
#include <vector>
#include <set>
#include <map>

#define MEDLOADERFORPV_EXPORTS
#include "MEDLoaderForPV.h"

class MEDLOADERFORPV_EXPORT ElectromagnetismRotationStatus
{
public:
  ElectromagnetismRotationStatus():_status(false) { }
  ElectromagnetismRotationStatus(const char *name);
  bool getStatus() const { return _status; }
  void setStatus(bool status) const { _status=status; }
  void cpyStatusFrom(const ElectromagnetismRotationStatus& other) { _status=other._status; }
  std::string getName() const { return _name; }
  void resetStatus() const { _status=false; }
  const char *getKeyOfEntry() const { return _ze_key_name.c_str(); }
  virtual void printMySelf(std::ostream& os) const;
  virtual bool isSameAs(const ElectromagnetismRotationStatus& other) const;
protected:
  mutable bool _status;
  std::string _name;
  std::string _ze_key_name;
};

class MEDLOADERFORPV_EXPORT ElectromagnetismRotationGrp : public ElectromagnetismRotationStatus
{
public:
  ElectromagnetismRotationGrp(const char *name);
  void setFamilies(const std::vector<std::string>& fams) { _fams=fams; }
  const std::vector<std::string>& getFamiliesLyingOn() const { return _fams; }
  bool isSameAs(const ElectromagnetismRotationGrp& other) const;
  static const char* start() { return "GRP_"; }
public:
  std::vector<std::string> _fams;
};

class MEDLOADERFORPV_EXPORT ElectromagnetismRotationFam : public ElectromagnetismRotationStatus
{
public:
  ElectromagnetismRotationFam(const char *name);
  void printMySelf(std::ostream& os) const;
  void fillIdsToKeep(std::set<int>& s) const;
  int getId() const { return _id; }
  bool isSameAs(const ElectromagnetismRotationFam& other) const;
  static const char* start() { return "FAM_"; }
private:
  int _id;
};

class vtkInformationDataObjectMetaDataKey;
class vtkMutableDirectedGraph;
class vtkInformation;

class MEDLOADERFORPV_EXPORT ElectromagnetismRotationInternal
{
public:
  void loadFrom(vtkMutableDirectedGraph *sil);
  int getNumberOfEntries() const;
  const char *getMeshName() const;
  const char *getKeyOfEntry(int i) const;
  bool getStatusOfEntryStr(const char *entry) const;
  void setStatusOfEntryStr(const char *entry, bool status);
  void printMySelf(std::ostream& os) const;
  std::vector<int> getFamiliesIdsOnGroup(const std::string& groupName) const;
  std::set<int> getIdsToKeep() const;
  std::vector< std::pair<std::string,std::vector<int> > > getAllGroups() const;
  void clearSelection() const;
  int getIdOfFamily(const std::string& famName) const;
  static bool IndependantIsInformationOK(vtkInformationDataObjectMetaDataKey *medReaderMetaData, vtkInformation *info);
private:
  std::map<std::string,int> computeFamStrIdMap() const;
  const ElectromagnetismRotationStatus& getEntry(const char *entry) const;
  ElectromagnetismRotationStatus& getEntry(const char *entry);
private:
  std::vector<ElectromagnetismRotationGrp> _groups;
  std::vector<ElectromagnetismRotationFam> _fams;
  mutable std::vector< std::pair<std::string,bool> > _selection;
  std::string _mesh_name;
};

