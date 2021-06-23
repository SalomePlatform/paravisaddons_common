from paraview.simple import *
from vtk.util import numpy_support
import numpy as np

def MyAssert(b):
    if not b:
        raise RuntimeError("Oooops : Assertion failed !")
NB_OF_LINES_FOR_1TS = 65 # not 66 because AB is not for inst == 2.0 !?
cr = ContactReader(registrationName='TEST.rco', FileName='23711_without_time.rco')
cr.UpdatePipelineInformation()
anim = GetAnimationScene()
tk = anim.TimeKeeper
times = tk.TimestepValues
MyAssert(times==[])
# time step 0
cr.UpdatePipeline()
ds0 = servermanager.Fetch(cr)
MyAssert([ds0.GetPointData().GetArrayName(i) for i in range(ds0.GetPointData().GetNumberOfArrays())]==["Resultante"])
MyAssert(ds0.GetNumberOfPoints()%NB_OF_LINES_FOR_1TS == 0)
step = int(ds0.GetNumberOfPoints()/NB_OF_LINES_FOR_1TS)
ts0 = numpy_support.vtk_to_numpy(ds0.GetPointData().GetAbstractArray(0))[::step]
MyAssert(ts0.shape == (NB_OF_LINES_FOR_1TS,3))
