# Copyright (C) 2021-2023  CEA, EDF
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

from medcoupling import *

data_file="Res_sisy.med"
poly_file="contour.med"
cutoff=0.

polygon=ReadMeshFromFile(poly_file)
data=MEDFileMesh.New(data_file)
mesh=data[0]
mesh.changeSpaceDimension(2,0.)
mesh=mesh.deepCopy()
assert(polygon.getCoords()[:,2].isUniform(0,1e-12))
polygon.changeSpaceDimension(2,0.)
polygon.mergeNodes(1e-12) #
polygon_1sgt=MEDCoupling1SGTUMesh(polygon)
conn=polygon_1sgt.getNodalConnectivity().deepCopy()
conn.rearrange(2)
notNullCellsPolygon=(conn[:,0]-conn[:,1]).findIdsNotEqual(0)
polygon=polygon_1sgt.buildUnstructured()[notNullCellsPolygon]
polygon_2d=MEDCouplingUMesh("mesh",2)
polygon_2d.setCoords(polygon.getCoords().deepCopy())
polygon_2d.allocateCells()
conn=MEDCoupling1SGTUMesh(polygon).getNodalConnectivity()
conn.rearrange(2)
conn2=conn.fromLinkedListOfPairToList()
assert(conn2[0]==conn2[-1])
conn2.popBackSilent()
polygon_2d.insertNextCell(NORM_POLYGON,conn2.getValues())
clockWise = polygon_2d.getMeasureField(False).getIJ(0,0) > 0.
#
side={True : 1 , False : 0}[clockWise]

(Xmin,Xmax),(Ymin,Ymax)=mesh.getBoundingBox()
TmpCenter=( (Xmin+Xmax)/2., (Ymin+Ymax)/2. )
Tmpalpha=1/max(Xmax-Xmin,Ymax-Ymin)
mesh.translate(-DataArrayDouble(TmpCenter,1,2))
mesh.scale([0.,0.],Tmpalpha)
polygon.translate(-DataArrayDouble(TmpCenter,1,2))
polygon.scale([0.,0.],Tmpalpha)
mesh2,line_inter,cellid_in_2d,cellid_in1d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(mesh,polygon,1e-12)
coo=mesh2.getCoords().deepCopy()
coo2=coo*(1/Tmpalpha)+TmpCenter
mesh2.setCoords(coo2)

twodcells_to_remove = cellid_in1d[:,side]
twodcells_to_keep = cellid_in1d[:,(side+1)%2]
ids=twodcells_to_keep.findIdsNotEqual(-1)
twodcells_to_keep=twodcells_to_keep[ids]
hotspot = twodcells_to_keep[0]

twodcells_to_keep = twodcells_to_keep[ids] # les cells2D de bord du domaine dans le referentiel de sortie 2D
twodcells_to_remove.sort()
ids=twodcells_to_remove.findIdsNotEqual(-1)
twodcells_to_remove=twodcells_to_remove[ids]
allcells=twodcells_to_remove.buildComplement(mesh2.getNumberOfCells())
mesh2_without_cells_around_polygon=mesh2[allcells]
grps=mesh2_without_cells_around_polygon.partitionBySpreadZone()
grps=[allcells[elt] for elt in grps]
assert(len(grps)==2)
zeGrp = None
if (hotspot in grps[0]) and (hotspot not in grps[1]):
    zeGrp = grps[0]
if (hotspot not in grps[0]) and (hotspot in grps[1]):
    zeGrp = grps[1]
if not zeGrp:
    raise RuntimeError("Ooops")
    pass
mesh3 = mesh2[zeGrp]
totVol = mesh3.getMeasureField(True).accumulate()[0]
refVol = polygon_2d.getMeasureField(True).accumulate()[0]
assert(abs((totVol-refVol)/refVol)<1e-6)
#
original_cell_ids_2d=cellid_in_2d[zeGrp] ; original_cell_ids_2d.sort() # les cells dans le referentiel original 2D
all_cut_2d_cells=cellid_in1d[:]
all_cut_2d_cells.rearrange(1)
all_cut_2d_cells.sort()
all_cut_2d_cells=all_cut_2d_cells.buildUnique() # les cells qui ont subit un split dans le referentiel de sortie 2D

all_cut_2d_cells_origin=cellid_in_2d[all_cut_2d_cells]
untouched_2d_cells=original_cell_ids_2d.buildSubstraction(all_cut_2d_cells_origin)

cells_at_boundary=all_cut_2d_cells.buildIntersection(zeGrp) # les cellules decoupées dans le referentiel de sortie 2D
cells_at_boundary_origin=cellid_in_2d[cells_at_boundary]
mesh3=mesh2[cells_at_boundary]
vol = mesh3.getMeasureField(True).getArray()
#volRef = mesh[cells_at_boundary_origin].getMeasureField(True).getArray()
#centers=mesh3.computeCellCenterOfMass()
#
tmp0=mesh[untouched_2d_cells]
(Xmin,Xmax),(Ymin,Ymax)=mesh3.getBoundingBox()
ZeCenter=( (Xmin+Xmax)/2., (Ymin+Ymax)/2. )
alpha=1/max(Xmax-Xmin,Ymax-Ymin)
mesh3.translate(-DataArrayDouble(ZeCenter,1,2))
mesh3.scale([0.,0.],alpha)
centers_around_zero=mesh3.computeCellCenterOfMass()
centers=centers_around_zero*(1/alpha)+ZeCenter
# 
evolution_multiTS=MEDFileFloatFieldMultiTS(data_file,"EVOLUTION",False)
for ts in range(10):
    evolution_1TS=evolution_multiTS[ts]
    evolution_1TS.loadArrays()
    f=evolution_1TS.field(data).convertToDblField()
    f_easy=f[untouched_2d_cells]
    f_easy.write("f_easy.med") #
    weights = f_easy.getDiscretization().getMeasureField(f_easy.getMesh(),True).getArray()
    positive_ids = f_easy.getArray().findIdsGreaterThan(cutoff)
    negative_ids = positive_ids.buildComplement(f_easy.getMesh().getNumberOfNodes())
    positive_part = (f_easy.getArray()[positive_ids]*weights[positive_ids]).accumulate(0)
    negative_part = (f_easy.getArray()[negative_ids]*weights[negative_ids]).accumulate(0)
    #
    f_hard = f[cells_at_boundary_origin]
    hard_part = f_hard.getValueOnMulti(centers)
    positive_ids_hard = hard_part.findIdsGreaterThan(cutoff)
    negative_ids_hard = positive_ids_hard.buildComplement(len(cells_at_boundary_origin))
    positive_part_hard=(hard_part[positive_ids_hard]*vol[positive_ids_hard]).accumulate(0)
    negative_part_hard=(hard_part[negative_ids_hard]*vol[negative_ids_hard]).accumulate(0)
    print("Time step %d (%ld)"%(ts,evolution_1TS.getTime()[-1]),positive_part+positive_part_hard+negative_part+negative_part_hard)
    evolution_1TS.unloadArrays()
    pass
