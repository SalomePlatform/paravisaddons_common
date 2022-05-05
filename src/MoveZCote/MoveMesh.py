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

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'MED Reader'
#f3d_gouttedomed = MEDReader(FileName='/home/H87074/TMP52/f3d_gouttedo.med')
#f3d_gouttedomed.AllArrays = ['TS0/MESH/ComSup0/COTE Z@@][@@P1', 'TS0/MESH/ComSup0/VITESSE U@@][@@P1', 'TS0/MESH/ComSup0/VITESSE V@@][@@P1', 'TS0/MESH/ComSup0/VITESSE W@@][@@P1']
#f3d_gouttedomed.AllTimeSteps = ['0000', '0001', '0002', '0003', '0004', '0005', '0006', '0007', '0008', '0009', '00010']

source = GetActiveSource()
renderView1 = GetActiveViewOrCreate('RenderView')
# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Calculator'
calculator1 = Calculator(Input=source)

# Properties modified on calculator1
calculator1.ResultArrayName = 'DisplacementsZ'
calculator1.Function = 'COTE Z-coordsZ'

# get color transfer function/color map for 'DisplacementsZ'
displacementsZLUT = GetColorTransferFunction('DisplacementsZ')

# show data in view
#calculator1Display = Show(calculator1, renderView1)

# hide data in view
Hide(source, renderView1)

# show color bar/color legend
#calculator1Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'DisplacementsZ'

# create a new 'Warp By Scalar'
warpByScalar1 = WarpByScalar(Input=calculator1)
warpByScalar1.Scalars = ['POINTS', 'DisplacementsZ']

