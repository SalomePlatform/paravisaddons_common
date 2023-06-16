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

from MEDLoader import *

fname="hydrau_test1.med"
meshName="mesh"
arr=DataArrayDouble([0,1,2,3,4,5])
m=MEDCouplingCMesh()
m.setCoords(arr,arr)
m=m.buildUnstructured()
m.setName(meshName)
m.simplexize(0)
WriteMesh(fname,m,True)
#
f=MEDCouplingFieldDouble(ON_NODES)
f.setMesh(m)
f.setName("Field")
arr=m.getCoords().magnitude()
f.setArray(arr)
for i in range(10):
    arr+=0.1
    f.setTime(float(i),i,0)
    WriteFieldUsingAlreadyWrittenMesh(fname,f)
    pass

