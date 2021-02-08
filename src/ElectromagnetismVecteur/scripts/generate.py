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

arr = DataArrayDouble([0,1,2])
m = MEDCouplingCMesh()
m.setCoords(arr,arr)
m = m.buildUnstructured()
#m.simplexize(0)
m.changeSpaceDimension(3,0.)
m.setName("mesh")
f = MEDCouplingFieldDouble(ON_CELLS)
f.setMesh(m)
f.setName("field")
f.setArray( DataArrayDouble([(1,1,0), (2,-2,0), (-3,3,0), (-4,-4,0)]) )
f.getArray().setInfoOnComponents(["X","Y","Z"])
f.checkConsistencyLight()
f.write("test.med")
f2 = f.deepCopy()
f2.setArray( DataArrayDouble([(1,0,0), (0,-2,0), (-3,0,0), (0,4,0)]) )
f2.setName("field2")
WriteFieldUsingAlreadyWrittenMesh("test.med",f2)
