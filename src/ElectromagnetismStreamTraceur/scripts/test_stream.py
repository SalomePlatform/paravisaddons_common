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

from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

def MyAssert(clue):
    if not clue:
        raise RuntimeError("Assertion failed !")

testmed = MEDReader(FileName='test.med')
testmed.AllArrays = ['TS0/mesh/ComSup0/field@@][@@P0', 'TS0/mesh/ComSup0/field2@@][@@P0', 'TS0/mesh/ComSup0/mesh@@][@@P0']
testmed.AllTimeSteps = ['0000']
streamTraceur1 = LigneDeChamp(Input=testmed,SeedType='Point Source')
streamTraceur1.SeedType.Radius = 1
streamTraceur1.SeedType.Center = [ 7.23,7.26,3.42 ]
streamTraceur1.Vectors = ['CELLS', "field"] #OrientationArray
streamTraceur1.UpdatePipeline()
ds0 = servermanager.Fetch(streamTraceur1)
MyAssert(ds0.GetNumberOfCells()==200)
