# Copyright (C) 2021  CEA/DEN, EDF R&D
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
import math

def ReturnInertia(p,OM,area_field):
    base_X = DataArrayDouble(len(OM),3) ; base_X[:]=p
    dist_to_base_X = (OM-DataArrayDouble.Dot(OM,base_X)*base_X).magnitude()
    inertia = (dist_to_base_X*dist_to_base_X*area_field).accumulate()[0]
    return inertia

def fffff(initialVect,normalFace,OM,area_field, posToIterate):
    li=[]
    for zePos in posToIterate:
        p = initialVect.deepCopy()
        MEDCouplingPointSet.Rotate3DAlg([0,0,0],normalFace.getValues(),zePos/float(180)*math.pi,p)
        inertia = ReturnInertia(p,OM,area_field)
        li.append((zePos,p.deepCopy(),inertia))
    return max(li,key=lambda x:x[2])

def fff(initialVect,normalFace,OM,area_field):
    pos = fffff(initialVect,normalFace,OM,area_field,[i*float(10) for i in range(18)])[0]
    for expo in range(5):
        pos,p,v = fffff(initialVect,normalFace,OM,area_field,[pos+i*(10**-expo) for i in range(-9,10)])
    return pos,p,v

"""
EDF24091 : probleme de calcul dans TorseurCIH
"""

fname = "Case_22_09_21.med"
fieldName = "LIN_____SIEF_NOEU"
mm=MEDFileMesh.New(fname)
m=mm[0]
f1ts = MEDFileField1TS(fname,fieldName)
f = f1ts.field(mm)
m = f.getMesh()
area_field = m.getMeasureField(True)
area = area_field.accumulate()[0] # 1
centerOfMassField = m.computeCellCenterOfMass()
centerOfMass = DataArrayDouble([elt/area for elt in (centerOfMassField*area_field.getArray()).accumulate()],1,3) # 2
m.unPolyze()
tri = MEDCoupling1SGTUMesh(m)
assert(tri.getCellModelEnum()==NORM_TRI3)
#
fCell = f.nodeToCellDiscretization()
(fCell.getArray()[:,[0,1,2]]*area).accumulate()
ids = area_field.getArray().findIdsLowerThan(1e-7).buildComplement(m.getNumberOfCells())
#fCell[ids]
eqn = m[ids].computePlaneEquationOf3DFaces()[:,:3]
eqn /= eqn.magnitude()
area_vector = eqn*area_field.getArray()[ids]
matrix = fCell[ids].getArray()
#
F_x = matrix[:,0]*eqn[:,0] + matrix[:,3]*eqn[:,1] + matrix[:,4]*eqn[:,2]
F_y = matrix[:,3]*eqn[:,0] + matrix[:,1]*eqn[:,1] + matrix[:,5]*eqn[:,2]
F_z = matrix[:,4]*eqn[:,0] + matrix[:,5]*eqn[:,1] + matrix[:,2]*eqn[:,2]
#
F = DataArrayDouble.Meld([F_x,F_y,F_z])
F[:] *= area_vector
#
ZeForce = DataArrayDouble(F.accumulate(),1,3)
normalFace = DataArrayDouble(eqn.accumulate(),1,3)
normalFace /= normalFace.magnitude()[0]
ForceNormale = DataArrayDouble.Dot(ZeForce,normalFace)[0]*normalFace # 3
TangentForce = ZeForce-ForceNormale # 4
#
bary = m[ids].computeCellCenterOfMass()
OM = bary-centerOfMass
momentum = DataArrayDouble(DataArrayDouble.CrossProduct(OM,F).accumulate(),1,3) # 5
# Inertie
InertiaNormale = (DataArrayDouble.Dot(OM,OM)*area_field.getArray()[ids]).accumulate()[0] # 6_A
base = DataArrayDouble(DataArrayDouble.GiveBaseForPlane(normalFace),3,3)
angle, tangentPrincipal, inertiaPrincipal = fff(base[0],normalFace,OM,area_field.getArray()[ids])
tangentOther = DataArrayDouble.CrossProduct(normalFace,tangentPrincipal)
inertiaOther = ReturnInertia(tangentOther,OM,area_field.getArray()[ids])
"""
base_X = DataArrayDouble(len(ids),3) ; base_X[:]=base[0]
base_Y = DataArrayDouble(len(ids),3) ; base_Y[:]=base[1]
dist_to_base_X = (OM-DataArrayDouble.Dot(OM,base_X)*base_X).magnitude()
dist_to_base_Y = (OM-DataArrayDouble.Dot(OM,base_Y)*base_Y).magnitude()
inertia_mat_0 = (dist_to_base_Y*dist_to_base_Y*area_field.getArray()[ids]).accumulate()[0]
inertia_mat_1 = (dist_to_base_X*dist_to_base_X*area_field.getArray()[ids]).accumulate()[0]
inertia_mat_01 = -(dist_to_base_X*dist_to_base_Y*area_field.getArray()[ids]).accumulate()[0]
from numpy import linalg as LA
import numpy as np
mat = np.matrix([[inertia_mat_0, inertia_mat_01], [inertia_mat_01, inertia_mat_1]])
v,w = LA.eig(mat)
pos_of_max = max([(i,elt) for (i,elt) in enumerate(v)],key=lambda x: x[1])[0]
u0 = DataArrayDouble(np.array(w[:,pos_of_max])) ; u0.rearrange(2)
v0 = DataArrayDouble(np.array(w[:,1-pos_of_max])) ; v0.rearrange(2)
#
I_new_base_0 = v[pos_of_max] # 6_B
new_base_0 = u0[0,0]*base[0]+u0[0,1]*base[1] # 6_B
#new_base_1 = v0[0,0]*base[0]+v0[0,1]*base[1]
new_base_1 = DataArrayDouble.CrossProduct(normalFace,new_base_0)
new_base_Y = DataArrayDouble(len(ids),3) ; new_base_Y[:]=new_base_1
new_dist_to_base_Y = (OM-DataArrayDouble.Dot(OM,new_base_Y)*new_base_Y).magnitude()
I_new_base_1 = (new_dist_to_base_Y*new_dist_to_base_Y*area_field.getArray()[ids]).accumulate()[0]
"""
"""
new_base_X = DataArrayDouble(len(ids),3) ; new_base_X[:]=new_base_0
new_dist_to_base_X = (OM-DataArrayDouble.Dot(OM,new_base_X)*new_base_X).magnitude()
I_new_base_0 = (new_dist_to_base_X*new_dist_to_base_X*area_field.getArray()[ids]).accumulate()[0]
tmp=m[ids] ; tmp.zipCoords()
f=MEDCouplingFieldDouble(ON_CELLS)
f.setMesh(tmp)
f.setArray(new_dist_to_base_X)
f.setName("dist")
f.writeVTK("distEig.vtu")"""
#mat*w[:,0]

print("{:g}".format(ForceNormale[0,0]))
