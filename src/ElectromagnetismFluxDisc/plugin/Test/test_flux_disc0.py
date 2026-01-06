# Copyright (C) 2021-2026  CEA, EDF
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

#  -*- coding: utf-8 -*-

from paraview.simple import *
from vtk.util import numpy_support
import numpy as np
import vtk

from medcoupling import *

def MyAssert(clue):
    if not clue:
        raise RuntimeError("Assertion failed !")

def Write(ds,fn):
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetInputData(ds)
    writer.SetFileName(fn)
    writer.Update()

def MCField(fn,fieldName):
    mm = MEDFileMesh.New(fn)
    fs = MEDFileFields.New(fn)
    return fs[fieldName][0].field(mm)

def ComputeFlux(mc_f):
    area = mc_f.getMesh().getMeasureField(True).getArray()
    mc_vector = DataArrayDouble(vector,1,3)
    mc_vector /= mc_vector.magnitude()[0]
    mc_vector = DataArrayDouble.Aggregate(mc_f.getNumberOfTuples()*[mc_vector])
    flux_per_cell = DataArrayDouble.Dot(mc_f.getArray(),mc_vector)*area
    return flux_per_cell.accumulate()[0]

ref_flux_b_a = -60.47
ref_flux_h_a = -127446
ref_flux_hsomega = 1896.991081794051
vector = [0.13,-0.47,0.87]
center = [0.0, 0.0, 3900.0]

mr = MEDReader(FileNames='/home/H87074/TMP81/paravistests_new/FilesForTests/Electromagnetism/mesh_benjamin_8_sept_2020.med')
mr.FieldsStatus = ['TS0/Mesh_1/ComSup0/B_A@@][@@P0', 'TS0/Mesh_1/ComSup0/B_HsOmega@@][@@P0', 'TS0/Mesh_1/ComSup0/H_A@@][@@P0', 'TS0/Mesh_1/ComSup0/H_HsOmega@@][@@P0']
mr.TimesFlagsStatus = ['0000', '0001', '0002', '0003', '0004', '0005', '0006', '0007', '0008', '0009']
mr.UpdatePipeline()

resRadial = 200
resCircum = 200
fluxDisc1 = FluxDisc(Input=mr)
fluxDisc1.SliceType = 'Cylinder'
fluxDisc1.SliceType.Center = center
fluxDisc1.SliceType.Axis = vector
fluxDisc1.SliceType.Radius = 293
fluxDisc1.RadialResolution = resRadial
fluxDisc1.CircumferentialResolution = resCircum

mw = MEDWriter(FileName="Base.med",Input=fluxDisc1)
mw.UpdatePipeline()

ds0 = servermanager.Fetch(fluxDisc1)
MyAssert(ds0.GetNumberOfCells()==resRadial*resCircum)
MyAssert(ds0.GetNumberOfPoints()==resRadial*resCircum+1)

flux_b_a = ds0.GetBlock(0).GetCellData().GetArray("B_A_flux").GetValue(0)
MyAssert(np.isclose(ref_flux_b_a,flux_b_a,0.1,0))

flux_h_a = ds0.GetBlock(0).GetCellData().GetArray("H_A_flux").GetValue(0)
#MyAssert(np.isclose(ref_flux_h_a,flux_h_a,0.01,0))

flux_hsomega = ds0.GetBlock(0).GetCellData().GetArray("H_HsOmega_flux").GetValue(0)
#MyAssert(np.isclose(ref_flux_hsomega,flux_hsomega,0.01,0))

Write(ds0.GetBlock(0),"Base.vtu")
mc_f_base = MCField("Base.med","H_HsOmega")

# On inverse l'axe et on verifie que le flux est inversé

fluxDisc1.SliceType.Axis = [-elt for elt in vector]
mw = MEDWriter(FileName="Invert.med",Input=fluxDisc1)
mw.UpdatePipeline()

ds1 = servermanager.Fetch(fluxDisc1)

MyAssert(ds1.GetNumberOfCells()==resRadial*resCircum)
MyAssert(ds1.GetNumberOfPoints()==resRadial*resCircum+1)

flux_b_a = ds1.GetBlock(0).GetCellData().GetArray("B_A_flux").GetValue(0)
#MyAssert(np.isclose(-ref_flux_b_a,flux_b_a,0.1,0))

flux_h_a = ds1.GetBlock(0).GetCellData().GetArray("H_A_flux").GetValue(0)
#MyAssert(np.isclose(-ref_flux_h_a,flux_h_a,0.1,0))

flux_hsomega = ds1.GetBlock(0).GetCellData().GetArray("H_HsOmega_flux").GetValue(0)
#MyAssert(np.isclose(-ref_flux_hsomega,flux_hsomega,0.1,0))


Write(ds0.GetBlock(0),"Invert.vtu")
mc_f_invert = MCField("Invert.med","H_HsOmega")

###

ComputeFlux(mc_f_base)
ComputeFlux(mc_f_invert)

# debug : pourquoi une telle difference de flux sur H_HsOmega ?

a  = mc_f_base.getMesh().getMeasureField(True).getArray()
b  = mc_f_invert.getMesh().getMeasureField(True).getArray()
assert(a.isEqual(b,1e-10))

# OK pour l'area

m_base = mc_f_base.getMesh()
m_invert = mc_f_invert.getMesh()
base_base = DataArrayDouble.GiveBaseForPlane(vector)

newX = DataArrayDouble(resRadial*resCircum+1,3) ; newX[:] = base_base[0]
newY = DataArrayDouble(resRadial*resCircum+1,3) ; newY[:] = base_base[1]
newZ = DataArrayDouble(resRadial*resCircum+1,3) ; newZ[:] = 0

base_new_coords=DataArrayDouble.Meld([DataArrayDouble.Dot(m_base.getCoords(),newX),DataArrayDouble.Dot(m_base.getCoords(),newY),newZ])
invert_new_coords=DataArrayDouble.Meld([DataArrayDouble.Dot(m_invert.getCoords(),newX),DataArrayDouble.Dot(m_invert.getCoords(),newY),newZ])

m_base.setCoords(base_new_coords)
m_invert.setCoords(invert_new_coords)

mc_f_base.write("base_dbg.vtu")
mc_f_invert.write("invert_dbg.vtu")

m_base.changeSpaceDimension(2,0.)
m_invert.changeSpaceDimension(2,0.)
a,b= m_base.getCellsContainingPoints(m_invert.computeCellCenterOfMass(),1e-12)
assert(b.isIota(len(b)))

a2 = a[:] ; a2.sort() ; assert( a2.isIota(len(a2) ) )

array_base_in_invert_ref = mc_f_base.getArray()[a]

not_of_ids = (array_base_in_invert_ref-mc_f_invert.getArray()).magnitude().findIdsGreaterThan(1.)
ok_ids = (array_base_in_invert_ref-mc_f_invert.getArray()).magnitude().findIdsLowerThan(1.)

ze_ids = ok_ids
ze_ids = not_of_ids
print(ComputeFlux(mc_f_invert[ze_ids]))
ze_ids_base_ref = a.findIdForEach(ze_ids)
print(ComputeFlux(mc_f_base[ze_ids_base_ref]))

mc_f_invert[ze_ids].write("invert_dbg_pb.vtu")
