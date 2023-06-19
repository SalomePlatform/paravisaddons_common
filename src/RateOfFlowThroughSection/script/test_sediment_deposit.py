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

# trace generated using paraview version 5.6.0-RC1-3-g7bafc2b

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'SerafinReader'
res_sisy_total_Tronqueres = SerafinReader(FileName='/home/H87074/TMP168_HYDRAU/17372_SEDIMENTS/Res_sisy_total_Tronque.res')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1268, 607]

# show data in view
res_sisy_total_TronqueresDisplay = Show(res_sisy_total_Tronqueres, renderView1)

# trace defaults for the display properties.
res_sisy_total_TronqueresDisplay.Representation = 'Surface'
res_sisy_total_TronqueresDisplay.ColorArrayName = [None, '']
res_sisy_total_TronqueresDisplay.OSPRayScaleArray = 'EVOLUTION      '
res_sisy_total_TronqueresDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
res_sisy_total_TronqueresDisplay.SelectOrientationVectors = 'EVOLUTION      '
res_sisy_total_TronqueresDisplay.ScaleFactor = 466.36875000000003
res_sisy_total_TronqueresDisplay.SelectScaleArray = 'EVOLUTION      '
res_sisy_total_TronqueresDisplay.GlyphType = 'Arrow'
res_sisy_total_TronqueresDisplay.GlyphTableIndexArray = 'EVOLUTION      '
res_sisy_total_TronqueresDisplay.GaussianRadius = 23.3184375
res_sisy_total_TronqueresDisplay.SetScaleArray = ['POINTS', 'EVOLUTION      ']
res_sisy_total_TronqueresDisplay.ScaleTransferFunction = 'PiecewiseFunction'
res_sisy_total_TronqueresDisplay.OpacityArray = ['POINTS', 'EVOLUTION      ']
res_sisy_total_TronqueresDisplay.OpacityTransferFunction = 'PiecewiseFunction'
res_sisy_total_TronqueresDisplay.DataAxesGrid = 'GridAxesRepresentation'
res_sisy_total_TronqueresDisplay.SelectionCellLabelFontFile = ''
res_sisy_total_TronqueresDisplay.SelectionPointLabelFontFile = ''
res_sisy_total_TronqueresDisplay.PolarAxes = 'PolarAxesRepresentation'
res_sisy_total_TronqueresDisplay.ScalarOpacityUnitDistance = 173.44832869720128

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
res_sisy_total_TronqueresDisplay.DataAxesGrid.XTitleFontFile = ''
res_sisy_total_TronqueresDisplay.DataAxesGrid.YTitleFontFile = ''
res_sisy_total_TronqueresDisplay.DataAxesGrid.ZTitleFontFile = ''
res_sisy_total_TronqueresDisplay.DataAxesGrid.XLabelFontFile = ''
res_sisy_total_TronqueresDisplay.DataAxesGrid.YLabelFontFile = ''
res_sisy_total_TronqueresDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
res_sisy_total_TronqueresDisplay.PolarAxes.PolarAxisTitleFontFile = ''
res_sisy_total_TronqueresDisplay.PolarAxes.PolarAxisLabelFontFile = ''
res_sisy_total_TronqueresDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
res_sisy_total_TronqueresDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'ShapeReader'
contourVolume1shp = ShapeReader(FileName='/home/H87074/TMP168_HYDRAU/17372_SEDIMENTS/ContourVolume1.shp')

# show data in view
contourVolume1shpDisplay = Show(contourVolume1shp, renderView1)

# trace defaults for the display properties.
contourVolume1shpDisplay.Representation = 'Surface'
contourVolume1shpDisplay.ColorArrayName = [None, '']
contourVolume1shpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
contourVolume1shpDisplay.SelectOrientationVectors = 'None'
contourVolume1shpDisplay.ScaleFactor = 28.764741869986757
contourVolume1shpDisplay.SelectScaleArray = 'None'
contourVolume1shpDisplay.GlyphType = 'Arrow'
contourVolume1shpDisplay.GlyphTableIndexArray = 'None'
contourVolume1shpDisplay.GaussianRadius = 1.4382370934993378
contourVolume1shpDisplay.SetScaleArray = [None, '']
contourVolume1shpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
contourVolume1shpDisplay.OpacityArray = [None, '']
contourVolume1shpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
contourVolume1shpDisplay.DataAxesGrid = 'GridAxesRepresentation'
contourVolume1shpDisplay.SelectionCellLabelFontFile = ''
contourVolume1shpDisplay.SelectionPointLabelFontFile = ''
contourVolume1shpDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
contourVolume1shpDisplay.DataAxesGrid.XTitleFontFile = ''
contourVolume1shpDisplay.DataAxesGrid.YTitleFontFile = ''
contourVolume1shpDisplay.DataAxesGrid.ZTitleFontFile = ''
contourVolume1shpDisplay.DataAxesGrid.XLabelFontFile = ''
contourVolume1shpDisplay.DataAxesGrid.YLabelFontFile = ''
contourVolume1shpDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contourVolume1shpDisplay.PolarAxes.PolarAxisTitleFontFile = ''
contourVolume1shpDisplay.PolarAxes.PolarAxisLabelFontFile = ''
contourVolume1shpDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
contourVolume1shpDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Transform'
transform1 = Transform(Input=contourVolume1shp)
transform1.Transform = 'Transform'

# Properties modified on transform1.Transform
transform1.Transform.Translate = [800000.0, 6500000.0, 0.0]

# Properties modified on transform1.Transform
transform1.Transform.Translate = [800000.0, 6500000.0, 0.0]

# show data in view
transform1Display = Show(transform1, renderView1)

# trace defaults for the display properties.
transform1Display.Representation = 'Surface'
transform1Display.ColorArrayName = [None, '']
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.SelectOrientationVectors = 'None'
transform1Display.ScaleFactor = 28.764741869986757
transform1Display.SelectScaleArray = 'None'
transform1Display.GlyphType = 'Arrow'
transform1Display.GlyphTableIndexArray = 'None'
transform1Display.GaussianRadius = 1.4382370934993378
transform1Display.SetScaleArray = [None, '']
transform1Display.ScaleTransferFunction = 'PiecewiseFunction'
transform1Display.OpacityArray = [None, '']
transform1Display.OpacityTransferFunction = 'PiecewiseFunction'
transform1Display.DataAxesGrid = 'GridAxesRepresentation'
transform1Display.SelectionCellLabelFontFile = ''
transform1Display.SelectionPointLabelFontFile = ''
transform1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
transform1Display.DataAxesGrid.XTitleFontFile = ''
transform1Display.DataAxesGrid.YTitleFontFile = ''
transform1Display.DataAxesGrid.ZTitleFontFile = ''
transform1Display.DataAxesGrid.XLabelFontFile = ''
transform1Display.DataAxesGrid.YLabelFontFile = ''
transform1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
transform1Display.PolarAxes.PolarAxisTitleFontFile = ''
transform1Display.PolarAxes.PolarAxisLabelFontFile = ''
transform1Display.PolarAxes.LastRadialAxisTextFontFile = ''
transform1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(contourVolume1shp, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(res_sisy_total_Tronqueres)

# create a new 'Sediment Deposit'
sedimentDeposit1 = SedimentDeposit(Input=res_sisy_total_Tronqueres,
    Source=transform1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [487587.4375, 6685895.75, 9562.2955669413]
renderView1.CameraFocalPoint = [487587.4375, 6685895.75, 0.0]
renderView1.CameraParallelScale = 2474.9042076238147

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
lineChartView1.ViewSize = [735, 607]
lineChartView1.ChartTitleFontFile = ''
lineChartView1.LeftAxisTitleFontFile = ''
lineChartView1.LeftAxisLabelFontFile = ''
lineChartView1.BottomAxisTitleFontFile = ''
lineChartView1.BottomAxisLabelFontFile = ''
lineChartView1.RightAxisLabelFontFile = ''
lineChartView1.TopAxisTitleFontFile = ''
lineChartView1.TopAxisLabelFontFile = ''

# get layout
layout1 = GetLayout()

# place view in the layout
layout1.AssignView(2, lineChartView1)

# show data in view
rateOfFlowThroughSection1Display = Show(sedimentDeposit1, lineChartView1)

# trace defaults for the display properties.
rateOfFlowThroughSection1Display.CompositeDataSetIndex = [0]
rateOfFlowThroughSection1Display.SeriesLabelPrefix = ''

# update the view to ensure updated data information
lineChartView1.Update()
