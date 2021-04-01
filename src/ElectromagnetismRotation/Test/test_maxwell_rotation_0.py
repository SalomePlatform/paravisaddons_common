from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
from vtk.util import numpy_support

def MyAssert(clue):
    if not clue:
        raise RuntimeError("Assertion failed !")

fileName="mesh_benjamin_8_sept_2020.med"
testmed = MEDReader(FileName=fileName)
testmed.AllArrays = ['TS0/Mesh_1/ComSup0/B_A@@][@@P0', 'TS0/Mesh_1/ComSup0/B_HsOmega@@][@@P0', 'TS0/Mesh_1/ComSup0/H_A@@][@@P0', 'TS0/Mesh_1/ComSup0/H_HsOmega@@][@@P0']
testmed.UpdatePipeline()

animationScene1 = GetAnimationScene()
timeKeeper1 = GetTimeKeeper()
animationScene1.UpdateAnimationUsingDataTimeSteps()

rog = RotationOfGroup(Input=testmed)
rog.AllGroups = ['GRP_airint_extruded', 'GRP_amortisseurs_extruded', 'GRP_arbre_extruded', 'GRP_inducteurGE_extruded', 'GRP_inducteurGS_extruded', 'GRP_inducteurPE_extruded', 'GRP_inducteurPS_extruded', 'GRP_mvt_extruded', 'GRP_rotor_extruded']
rog.AngularStep = "180+180"
rog.UpdatePipeline()
ds0 = servermanager.Fetch(rog)
MyAssert(ds0.GetBlock(0).GetNumberOfCells() == 3423)
rog.UpdatePipeline()
animationScene1.GoToNext()
rog.UpdatePipeline()
animationScene1.GoToNext()
rog.UpdatePipeline()
