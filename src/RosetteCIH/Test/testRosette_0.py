#  -*- coding: utf-8 -*-
# Copyright (C) 2025-2026  CEA, EDF
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

from paraview.simple import *
from pathlib import Path
import medcoupling as mc
import tempfile

arr = mc.DataArrayDouble(3) ; arr.iota()
m = mc.MEDCouplingCMesh() ; m.setCoords(arr,arr)
m = m.buildUnstructured()
m.changeSpaceDimension(3,0.)
m.setName("mesh")
f = mc.MEDCouplingFieldDouble(mc.ON_CELLS)
f.setName("SIRO_ELEM_T1")
f.setMesh(m)
arr=mc.DataArrayDouble([(1,0,0,1),(0,1,0,1),(-1,0,0,1),(0,-1,0,1)])
f.setArray(arr)
arr.setInfoOnComponents(["SIG_T1X","SIG_T1Y","SIG_T1Z","SIG_T1"])

with tempfile.TemporaryDirectory() as di:
    fname = "{}".format( Path(di) / "testRosette.med" )
    f.write(fname)

    rosette = MEDReader(registrationName='rosette.med', FileNames=[fname])
    rosette.FieldsStatus = ['TS0/mesh/ComSup0/SIRO_ELEM_T1@@][@@P0']

    rosettesdecontrainte1 = Rosettesdecontrainte(registrationName='Rosettesdecontrainte1', Input=rosette)
    rosettesdecontrainte1.TypeOfDisplay = "T1 only"
    rosettesdecontrainte1.ScaleFactor = 0.5
    rosettesdecontrainte1.WidthFactor = 0.1
    rosettesdecontrainte1.UpdatePipeline()

    res = servermanager.Fetch( rosettesdecontrainte1  )
    assert( res.GetBlock(0).GetNumberOfCells() == 8 )
