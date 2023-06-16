# Copyright (C) 2021-2023  CEA/DEN, EDF R&D
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

resu = "resu.med"
suffix_section = "5"
section = "section_{}.med".format(suffix_section)

print("Working with {}".format(section))
mm=MEDFileMesh.New(resu)
m=mm[0]
sec=ReadMeshFromFile(section)
assert(m.getCoords()[:,2].isUniform(0,1e-12))
m.setCoords(m.getCoords()[:,[0,1]])
assert(sec.getCoords()[:,2].isUniform(0,1e-12))
sec.setCoords(sec.getCoords()[:,[0,1]])
sec.mergeNodes(1e-12)
sec.zipCoords()
sec.removeDegenerated1DCells()
_,line_inter,cellid_in_2d,cellid_in1d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m,sec,1e-6)
line_inter.zipCoords()
#
TwoDcells=DataArrayInt(line_inter.getNumberOfCells())
for i,t in enumerate(cellid_in1d):
    candidates = [elt for elt in list(t) if elt != -1]
    if len(candidates)==0:
        TwoDcells[i]=-1
    TwoDcells[i]=candidates[0]
    pass
notFreeStyle1DCells = TwoDcells.findIdsNotEqual(-1)
n2oCells = TwoDcells[notFreeStyle1DCells]
TwoDcells=cellid_in_2d[n2oCells]
#
effective_line1d=line_inter[notFreeStyle1DCells]
# effective_2d_cells - maillage contenant pour chaque cellule du maillage 1D coupé la cellule 2D de resu qui la contient
effective_2d_cells=m[TwoDcells]
o2n=effective_2d_cells.zipCoordsTraducer()
n2o=o2n.invertArrayO2N2N2O(effective_2d_cells.getNumberOfNodes())
#
effective_line1d=MEDCoupling1SGTUMesh(effective_line1d) # change format of umesh to ease alg
effective_2d_cells=MEDCoupling1SGTUMesh(effective_2d_cells) # change format of umesh to ease alg
assert(effective_2d_cells.getCellModelEnum()==NORM_TRI3)
#
conn1d=effective_line1d.getNodalConnectivity()[:] ; conn1d.rearrange(2)
conn2d=effective_2d_cells.getNodalConnectivity()[:] ; conn2d.rearrange(3)
coo1d=effective_line1d.getCoords()
coo2d=effective_2d_cells.getCoords()
assert(len(conn2d)==len(conn1d))
# coeffs coeffs_integ and n2o are elements for matrix
h_water_mts=MEDFileFloatFieldMultiTS(resu,"HAUTEUR D\'EAU",False)
speed_mts=MEDFileFloatFieldMultiTS(resu,"VITESSE",False)
assert(len(h_water_mts.getPflsReallyUsed())==0)
assert(len(speed_mts.getPflsReallyUsed())==0)

h_out=DataArrayDouble(effective_line1d.getNumberOfCells()) ; h_out[:]=0.
v_out=DataArrayDouble(effective_line1d.getNumberOfCells()) ; v_out[:]=0.
# on calcule la matrice qui pour chaque cellule du la line 1D decoupee, donne 
# la contribution de chacun des nodes de la cell 2D a laquelle elle appartient.
matrix=effective_line1d.getNumberOfCells()*[None]

for i,(t1,t2) in enumerate(zip(conn1d,conn2d)):
    seg2=coo1d[list(t1)]
    tri3=coo2d[list(t2)]
    baryInfo,length=DataArrayDouble.ComputeIntegralOfSeg2IntoTri3(seg2,tri3)
    matrix[i]=[(n2o[i],j) for i,j in zip(list(t2),baryInfo)]
    pass
ortho=effective_line1d.buildUnstructured().buildOrthogonalField().getArray()

for ts in range(1):
    coeffs_integ=DataArrayDouble(effective_2d_cells.getNumberOfNodes()) ; coeffs_integ[:]=0
    h_water=h_water_mts[ts] ; h_water.loadArrays()
    speed=speed_mts[ts] ; speed.loadArrays()
    h_water_arr=h_water.getUndergroundDataArray()
    speed_arr=speed.getUndergroundDataArray().convertToDblArr()
    assert(speed_arr[:,2].isUniform(0,1e-12))
    speed_arr=speed_arr[:,[0,1]]
    for i in range(effective_line1d.getNumberOfCells()):
        row=matrix[i]
        h_out[i]=sum([b*h_water_arr[a] for a,b in row])
        speed=sum([b*speed_arr[a] for a,b in row])
        v_out[i] = float(DataArrayDouble.Dot(speed,ortho[i])[0])
        pass
    zeValue = abs((h_out*v_out*effective_line1d.getMeasureField(True).getArray()).accumulate()[0])
    print("ts %d (%d) = %lf"%(ts,int(h_water.getTime()[-1]),zeValue))
    pass
# Bug 21733
# Avant 1 == 1487.5
#       5 == 1434.8
# Apres 1 == 1600.814665
#       5 == 1600.837195
"""h_f=MEDCouplingFieldDouble(ON_CELLS) ; h_f.setMesh(effective_line1d)
h_f.setArray(h_out)
h_f.setName("HAUTEUR")
#
v_f=MEDCouplingFieldDouble(ON_CELLS) ; v_f.setMesh(effective_line1d)
v_f.setArray(v_out)
v_f.setName("VITESSE")
effective_line1d.write("line1d.med")
WriteFieldUsingAlreadyWrittenMesh("line1d.med",h_f)
WriteFieldUsingAlreadyWrittenMesh("line1d.med",v_f)"""

"""# avec calcul_2
ts 0 (0) = 1606.649455
ts 1 (7200) = 1534.516771
ts 2 (14400) = 1549.476531
ts 3 (21600) = 1551.205389
ts 4 (28800) = 1550.100327
ts 5 (36000) = 1547.519873
ts 6 (43200) = 1542.625840
ts 7 (50400) = 1540.418416
ts 8 (57600) = 1539.691491
ts 9 (64800) = 1542.502136
ts 10 (72000) = 1536.397618
ts 11 (79200) = 1536.609661
ts 12 (86400) = 1535.983922
ts 13 (93600) = 1537.728434
ts 14 (100800) = 1537.462885
ts 15 (108000) = 1537.290268
ts 16 (115200) = 1537.143315
ts 17 (122400) = 1537.037729
ts 18 (129600) = 1536.967132
ts 19 (136800) = 1536.924427
ts 20 (144000) = 1536.905037"""
