# Copyright (C) 2021-2022  CEA/DEN, EDF R&D
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

from MEDLoader import *
import numpy as np

fname="quadratic.med"
meshName="mesh"
coo=DataArrayDouble([(0,-1,0),(1,-1,0),(2,-1,0),#0->3
                     (-1,0,0),(0,0,0),(1,0,0),(2,0,0),#3->7
                     (2,1,0),(1,1,0),(0,1,0),#7->10
                     (-1,0,1),(0,0,1),(1,0,1),#10->13
                     (0,1,1),(1,1,1)#
                 ])
m0=MEDCouplingUMesh(meshName,3) ; m0.setCoords(coo)
m0.allocateCells()
m0.insertNextCell(NORM_TETRA4,[5,6,7,12])
m0.insertNextCell(NORM_PENTA6,[3,9,4,10,13,11])
m0.insertNextCell(NORM_HEXA8,[4,9,8,5,11,13,14,12]) ; m0.zipCoords()
m0.convertLinearCellsToQuadratic()
#
m1=MEDCouplingUMesh(meshName,2) ; m1.setCoords(coo)
m1.allocateCells()
m1.insertNextCell(NORM_TRI3,[3,4,0])
m1.insertNextCell(NORM_QUAD4,[0,4,5,1]) ; m1.zipCoords()
m1.convertLinearCellsToQuadratic()
#
m2=MEDCouplingUMesh(meshName,1) ; m2.setCoords(coo)
m2.allocateCells()
m2.insertNextCell(NORM_SEG2,[1,2]) ; m2.zipCoords()
m2.convertLinearCellsToQuadratic()
#
coos=[m0.getCoords(),m1.getCoords(),m2.getCoords()]
ncoos=[len(cooElt) for cooElt in coos]
arr=DataArrayDouble.Aggregate(coos)
a,b=arr.findCommonTuples(1e-12)
o2n,newNbNodes=DataArrayInt.ConvertIndexArrayToO2N(len(arr),a,b)
newArr=arr[o2n.invertArrayO2N2N2O(newNbNodes)]
#
m0.renumberNodesInConn(o2n) ; m0.setCoords(newArr)
m1.renumberNodesInConn(o2n[sum(ncoos[:1]):]) ; m1.setCoords(newArr)
m2.renumberNodesInConn(o2n[sum(ncoos[:2]):]) ; m2.setCoords(newArr)
a,b=newArr.areIncludedInMe(coo,1e-12)
assert(a) ; b.sort()
nodeIdsNotLinear=b.buildComplement(len(newArr))
#
mm=MEDFileUMesh()
mm[0]=m0 ; mm[-1]=m1 ; mm[-2]=m2
mm.write(fname,2)

#
centerOfCloud=[newArr[:,i].accumulate(0)/len(newArr) for i in range(newArr.getNumberOfComponents())]
farr=(newArr-centerOfCloud).magnitude()
farr[nodeIdsNotLinear]=0.
zeMax=farr.getMaxValue()[0]
#

fieldName="FieldDouble"
ff=MEDFileField1TS()
f=MEDCouplingFieldDouble(ON_NODES) ; f.setMesh(m0) ; f.setArray(farr) ; f.setName(fieldName) ; f.setTime(0.,0,0)
ff.setFieldNoProfileSBT(f)
ff.write(fname,0)

fieldName="FieldFloat"
ff=MEDFileFloatField1TS()
farr = DataArrayFloat(np.ones((len(farr)), dtype=np.float32))
f=MEDCouplingFieldFloat(ON_NODES) ; f.setMesh(m0) ; f.setArray(farr) ; f.setName(fieldName) ; f.setTime(0.,0,0)
ff.setFieldNoProfileSBT(f)
ff.write(fname,0)

fieldName="FieldInt32"
ff=MEDFileInt32Field1TS()
farr = DataArrayInt32(np.ones((len(farr)), dtype=np.int32))
f=MEDCouplingFieldInt32(ON_NODES) ; f.setMesh(m0) ; f.setArray(farr) ; f.setName(fieldName) ; f.setTime(0.,0,0)
ff.setFieldNoProfileSBT(f)
ff.write(fname,0)

fieldName="FieldInt64"
ff=MEDFileInt64Field1TS()
farr = DataArrayInt64(np.ones((len(farr)), dtype=np.int64))
f=MEDCouplingFieldInt64(ON_NODES) ; f.setMesh(m0) ; f.setArray(farr) ; f.setName(fieldName) ; f.setTime(0.,0,0)
ff.setFieldNoProfileSBT(f)
ff.write(fname,0)







#fieldName="FieldDouble"
#f=MEDCouplingFieldDouble(ON_NODES, ONE_TIME) ; f.setMesh(m0) ; f.setArray(farr) ; f.setName(fieldName) ; f.setTime(0.,0,0)
#farr-=zeMax ; farr.abs() ; farr[nodeIdsNotLinear]=0.
#farr.reverse()
#WriteFieldUsingAlreadyWrittenMesh(fname, f)
#
#fieldName="FieldFloat"
#farr = DataArrayFloat(np.ones((len(farr)), dtype=np.float32))
#f=MEDCouplingFieldFloat(ON_NODES, ONE_TIME) ; f.setMesh(m0) ; f.setArray(farr) ; f.setName(fieldName) ; f.setTime(0.,0,0)
#WriteFieldUsingAlreadyWrittenMesh(fname, f)
#
#fieldName="FieldInt32"
#farr = DataArrayInt32(np.ones((len(farr)), dtype=np.int32))
#f=MEDCouplingFieldInt32(ON_NODES, ONE_TIME) ; f.setMesh(m0) ; f.setArray(farr) ; f.setName(fieldName) ; f.setTime(0.,0,0)
#WriteFieldUsingAlreadyWrittenMesh(fname, f)

#fieldName="FieldInt64"
#farr = DataArrayInt64(np.ones((len(farr)), dtype=np.int64))
#f=MEDCouplingFieldInt64(ON_NODES, ONE_TIME) ; f.setMesh(m0) ; f.setArray(farr) ; f.setName(fieldName) ; f.setTime(0.,0,0)
#WriteField(fname, f, True)

from paraview.simple import *

quadraticmed = MEDReader(registrationName=fname, FileName=fname)
quadraticmed.AllArrays = ['TS0/mesh/ComSup0/FieldDouble@@][@@P1', 'TS0/mesh/ComSup0/FieldFloat@@][@@P1', 'TS0/mesh/ComSup0/FieldInt32@@][@@P1']
quadraticmed.AllTimeSteps = ['0000']

quadraticmed.UpdatePipeline()

quadraticToLinear1 = QuadraticToLinear(registrationName='QuadraticToLinear1', Input=quadraticmed)
quadraticToLinear1.UpdatePipeline()
