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
// Author : Anthony Geay (EDF R&D)

#include "ElectromagnetismRotationHelper.h"

#include "InterpKernelException.hxx"

#include "vtkInformation.h"
#include "vtkInformationDataObjectMetaDataKey.h"
#include "vtkAdjacentVertexIterator.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkDataSetAttributes.h"
#include "vtkStringArray.h"

#include <cstring>
#include <limits>
#include <algorithm>

const char ZE_SEPP[]="@@][@@";

ElectromagnetismRotationStatus::ElectromagnetismRotationStatus(const char *name):_status(false),_name(name)
{
}

void ElectromagnetismRotationStatus::printMySelf(std::ostream& os) const
{
  os << "      -" << _ze_key_name << "(";
  if(_status)
    os << "X";
  else
    os << " ";
  os << ")" << std::endl;
}

bool ElectromagnetismRotationStatus::isSameAs(const ElectromagnetismRotationStatus& other) const
{
  return _name==other._name && _ze_key_name==other._ze_key_name;
}

bool ElectromagnetismRotationGrp::isSameAs(const ElectromagnetismRotationGrp& other) const
{
  bool ret(ElectromagnetismRotationStatus::isSameAs(other));
  if(ret)
    return _fams==other._fams;
  else
    return false;
}

ElectromagnetismRotationGrp::ElectromagnetismRotationGrp(const char *name):ElectromagnetismRotationStatus(name)
{
  std::ostringstream oss; oss << ElectromagnetismRotationGrp::start() << name; _ze_key_name=oss.str();
}

ElectromagnetismRotationFam::ElectromagnetismRotationFam(const char *name):ElectromagnetismRotationStatus(name),_id(0)
{
  std::size_t pos(_name.find(ZE_SEPP));
  std::string name0(_name.substr(0,pos)),name1(_name.substr(pos+strlen(ZE_SEPP)));
  std::istringstream iss(name1);
  iss >> _id;
  std::ostringstream oss; oss << ElectromagnetismRotationFam::start() << name; _ze_key_name=oss.str(); _name=name0;
}

bool ElectromagnetismRotationFam::isSameAs(const ElectromagnetismRotationFam& other) const
{
  bool ret(ElectromagnetismRotationStatus::isSameAs(other));
  if(ret)
    return _id==other._id;
  else
    return false;
}

void ElectromagnetismRotationFam::printMySelf(std::ostream& os) const
{
  os << "      -" << _ze_key_name << " famName : \"" << _name << "\" id : " << _id << " (";
  if(_status)
    os << "X";
  else
    os << " ";
  os << ")" << std::endl;
}

void ElectromagnetismRotationFam::fillIdsToKeep(std::set<int>& s) const
{
  s.insert(_id);
}
///////////////////

bool ElectromagnetismRotationInternal::IndependantIsInformationOK(vtkInformationDataObjectMetaDataKey *medReaderMetaData, vtkInformation *info)
{
  // Check the information contain meta data key
  if(!info->Has(medReaderMetaData))
    return false;

  // Recover Meta Data
  vtkMutableDirectedGraph *sil(vtkMutableDirectedGraph::SafeDownCast(info->Get(medReaderMetaData)));
  if(!sil)
    return false;
  int idNames(0);
  vtkAbstractArray *verticesNames(sil->GetVertexData()->GetAbstractArray("Names",idNames));
  vtkStringArray *verticesNames2(vtkStringArray::SafeDownCast(verticesNames));
  if(!verticesNames2)
    return false;
  for(int i=0;i<verticesNames2->GetNumberOfValues();i++)
    {
      vtkStdString &st(verticesNames2->GetValue(i));
      if(st=="MeshesFamsGrps")
        return true;
    }
  return false;
}

const char *ElectromagnetismRotationInternal::getMeshName() const
{
  return this->_mesh_name.c_str();
}

void ElectromagnetismRotationInternal::loadFrom(vtkMutableDirectedGraph *sil)
{
  std::vector<ElectromagnetismRotationGrp> oldGrps(_groups); _groups.clear();
  std::vector<ElectromagnetismRotationFam> oldFams(_fams); _fams.clear();
  int idNames(0);
  vtkAbstractArray *verticesNames(sil->GetVertexData()->GetAbstractArray("Names",idNames));
  vtkStringArray *verticesNames2(vtkStringArray::SafeDownCast(verticesNames));
  vtkIdType id0;
  bool found(false);
  for(int i=0;i<verticesNames2->GetNumberOfValues();i++)
    {
      vtkStdString &st(verticesNames2->GetValue(i));
      if(st=="MeshesFamsGrps")
        {
          id0=i;
          found=true;
        }
    }
  if(!found)
    throw INTERP_KERNEL::Exception("There is an internal error ! The tree on server side has not the expected look !");
  vtkAdjacentVertexIterator *it0(vtkAdjacentVertexIterator::New());
  sil->GetAdjacentVertices(id0,it0);
  int kk(0),ll(0);
  while(it0->HasNext())
    {
      vtkIdType id1(it0->Next());
      std::string meshName(verticesNames2->GetValue(id1));
      this->_mesh_name=meshName;
      vtkAdjacentVertexIterator *it1(vtkAdjacentVertexIterator::New());
      sil->GetAdjacentVertices(id1,it1);
      vtkIdType idZeGrps(it1->Next());//zeGroups
      vtkAdjacentVertexIterator *itGrps(vtkAdjacentVertexIterator::New());
      sil->GetAdjacentVertices(idZeGrps,itGrps);
      while(itGrps->HasNext())
        {
          vtkIdType idg(itGrps->Next());
          ElectromagnetismRotationGrp grp(verticesNames2->GetValue(idg).c_str());
          vtkAdjacentVertexIterator *itGrps2(vtkAdjacentVertexIterator::New());
          sil->GetAdjacentVertices(idg,itGrps2);
          std::vector<std::string> famsOnGroup;
          while(itGrps2->HasNext())
            {
              vtkIdType idgf(itGrps2->Next());
              famsOnGroup.push_back(std::string(verticesNames2->GetValue(idgf)));
            }
          grp.setFamilies(famsOnGroup);
          itGrps2->Delete();
          _groups.push_back(grp);
        }
      itGrps->Delete();
      vtkIdType idZeFams(it1->Next());//zeFams
      it1->Delete();
      vtkAdjacentVertexIterator *itFams(vtkAdjacentVertexIterator::New());
      sil->GetAdjacentVertices(idZeFams,itFams);
      while(itFams->HasNext())
        {
          vtkIdType idf(itFams->Next());
          ElectromagnetismRotationFam fam(verticesNames2->GetValue(idf).c_str());
          _fams.push_back(fam);
        }
      itFams->Delete();
    }
  it0->Delete();
  // filter groups on cells
  std::vector<ElectromagnetismRotationGrp> groupsToKeep;
  std::size_t ii(0);
  for(auto grp : _groups)
  {
    std::vector<int> famIds(this->getFamiliesIdsOnGroup(grp.getName()));
    if ( std::all_of(famIds.begin(), famIds.end(), [](int i){ return i<0; }) )
      groupsToKeep.emplace_back(std::move(grp));
  }
  _groups = std::move(groupsToKeep);
  //
  std::size_t szg(_groups.size()),szf(_fams.size());
  if(szg==oldGrps.size() && szf==oldFams.size())
    {
      bool isSame(true);
      for(std::size_t i=0;i<szg && isSame;i++)
        isSame=_groups[i].isSameAs(oldGrps[i]);
      for(std::size_t i=0;i<szf && isSame;i++)
        isSame=_fams[i].isSameAs(oldFams[i]);
      if(isSame)
        {
          for(std::size_t i=0;i<szg;i++)
            _groups[i].cpyStatusFrom(oldGrps[i]);
          for(std::size_t i=0;i<szf;i++)
            _fams[i].cpyStatusFrom(oldFams[i]);
        }
    }
}

int ElectromagnetismRotationInternal::getNumberOfEntries() const
{
  std::size_t sz0(_groups.size());
  return (int)(sz0);
}

const char *ElectromagnetismRotationInternal::getKeyOfEntry(int i) const
{
  return _groups[i].getKeyOfEntry();
}

bool ElectromagnetismRotationInternal::getStatusOfEntryStr(const char *entry) const
{
  const ElectromagnetismRotationStatus& elt(getEntry(entry));
  return elt.getStatus();
}

void ElectromagnetismRotationInternal::setStatusOfEntryStr(const char *entry, bool status)
{
  _selection.emplace_back(entry,status);
}

const ElectromagnetismRotationStatus& ElectromagnetismRotationInternal::getEntry(const char *entry) const
{
  std::string entryCpp(entry);
  for(std::vector<ElectromagnetismRotationGrp>::const_iterator it0=_groups.begin();it0!=_groups.end();it0++)
    if(entryCpp==(*it0).getKeyOfEntry())
      return *it0;
  std::ostringstream oss; oss << "vtkElectromagnetismRotationInternal::getEntry : no such entry \"" << entry << "\"!";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

ElectromagnetismRotationStatus& ElectromagnetismRotationInternal::getEntry(const char *entry)
{
  std::string entryCpp(entry);
  for(std::vector<ElectromagnetismRotationGrp>::iterator it0=_groups.begin();it0!=_groups.end();it0++)
    if(entryCpp==(*it0).getKeyOfEntry())
      return *it0;
  std::ostringstream oss; oss << "vtkElectromagnetismRotationInternal::getEntry : no such entry \"" << entry << "\"!";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

void ElectromagnetismRotationInternal::printMySelf(std::ostream& os) const
{
  os << "Groups :" << std::endl;
  for(std::vector<ElectromagnetismRotationGrp>::const_iterator it0=_groups.begin();it0!=_groups.end();it0++)
    (*it0).printMySelf(os);
}

std::vector<int> ElectromagnetismRotationInternal::getFamiliesIdsOnGroup(const std::string& groupName) const
{
  for(auto grp : _groups)
  {
    if(grp.getName() == groupName)
    {
      std::vector<std::string> fams(grp.getFamiliesLyingOn());
      auto sz(fams.size());
      std::vector<int> famIds(sz);
      for(auto i = 0 ; i < sz ; ++i)
        famIds[i] = this->getIdOfFamily(fams[i]);
      return famIds;
    }
  }
  std::ostringstream oss; oss << "vtkElectromagnetismRotationInternal::getFamiliesIdsOnGroup : no such group \"" << groupName << "\"!";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

int ElectromagnetismRotationInternal::getIdOfFamily(const std::string& famName) const
{
  for(std::vector<ElectromagnetismRotationFam>::const_iterator it=_fams.begin();it!=_fams.end();it++)
    {
      if((*it).getName()==famName)
        return (*it).getId();
    }
  return std::numeric_limits<int>::max();
}

std::set<int> ElectromagnetismRotationInternal::getIdsToKeep() const
{
  for(auto it: _selection)
    {
      const ElectromagnetismRotationStatus& elt(getEntry(it.first.c_str()));
      elt.setStatus(it.second);
    }
  std::map<std::string,int> m(this->computeFamStrIdMap());
  std::set<int> s;
  for(std::vector<ElectromagnetismRotationGrp>::const_iterator it0=_groups.begin();it0!=_groups.end();it0++)
    {
      if((*it0).getStatus())
        {
          const std::vector<std::string>& fams((*it0).getFamiliesLyingOn());
          for(std::vector<std::string>::const_iterator it1=fams.begin();it1!=fams.end();it1++)
            {
              std::map<std::string,int>::iterator it2(m.find((*it1)));
              if(it2!=m.end())
                s.insert((*it2).second);
            }
        }
     }
  for(std::vector<ElectromagnetismRotationFam>::const_iterator it0=_fams.begin();it0!=_fams.end();it0++)
    if((*it0).getStatus())
      (*it0).fillIdsToKeep(s);
  return s;
}

// see reference : https://en.cppreference.com/w/cpp/iterator/iterator
class FamilyIterator : public std::iterator< std::input_iterator_tag, long, long, int*, int >
{
  long _num = 0;
  const ElectromagnetismRotationInternal *_egi = nullptr;
  const std::vector<std::string> *_fams = nullptr;
public:
  explicit FamilyIterator(long num , const ElectromagnetismRotationInternal *egi, const std::vector<std::string>& fams) : _num(num),_egi(egi),_fams(&fams) {}
  FamilyIterator& operator++() { ++_num; return *this;}
  bool operator==(const FamilyIterator& other) const {return _num == other._num;}
  bool operator!=(const FamilyIterator& other) const {return !(*this == other);}
  reference operator*() const {return _egi->getIdOfFamily((*_fams)[_num]);}
};

std::vector< std::pair<std::string,std::vector<int> > > ElectromagnetismRotationInternal::getAllGroups() const
{
    std::vector< std::pair<std::string,std::vector<int> > > ret;
    for(const auto&  grp : _groups)
    {
        const std::vector<std::string>& fams(grp.getFamiliesLyingOn());
        std::vector<int> famIds(FamilyIterator(0,this,fams),FamilyIterator(fams.size(),this,fams));
        if ( std::all_of(famIds.begin(), famIds.end(), [](int i){ return i<0; }) )// only groups on cells considered here
        {
          std::pair<std::string,std::vector<int> > elt(grp.getName(),std::move(famIds));
          ret.emplace_back(std::move(elt));
        }
    }
    return ret;
}

void ElectromagnetismRotationInternal::clearSelection() const
{
  _selection.clear();
  for(auto it : _groups)
    it.resetStatus();
  for(auto it : _fams)
    it.resetStatus();
}

std::map<std::string,int> ElectromagnetismRotationInternal::computeFamStrIdMap() const
{
  std::map<std::string,int> ret;
  for(std::vector<ElectromagnetismRotationFam>::const_iterator it0=_fams.begin();it0!=_fams.end();it0++)
    ret[(*it0).getName()]=(*it0).getId();
  return ret;
}
