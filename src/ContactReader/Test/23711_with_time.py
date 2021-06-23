from paraview.simple import *
from vtk.util import numpy_support
import numpy as np

def MyAssert(b):
    if not b:
        raise RuntimeError("Oooops : Assertion failed !")
NB_OF_LINES_FOR_1TS = 22
cr = ContactReader(registrationName='TEST.rco', FileName='23711_with_time.rco')
cr.UpdatePipelineInformation()
anim = GetAnimationScene()
tk = anim.TimeKeeper
times = tk.TimestepValues
MyAssert(times==[0.0,1.0,2.0])
Show(cr)
meshViewer = GetActiveView()
# time step 0
meshViewer.ViewTime = times[0]
meshViewer.Representations[0].UpdatePipeline()
Render()
cr.UpdatePipeline()
ds0 = servermanager.Fetch(cr)
MyAssert([ds0.GetPointData().GetArrayName(i) for i in range(ds0.GetPointData().GetNumberOfArrays())]==["Resultante"])
MyAssert(ds0.GetNumberOfPoints()%NB_OF_LINES_FOR_1TS == 0)
step = int(ds0.GetNumberOfPoints()/NB_OF_LINES_FOR_1TS)
ts0 = numpy_support.vtk_to_numpy(ds0.GetPointData().GetAbstractArray(0))[::step]
MyAssert(ts0.shape == (NB_OF_LINES_FOR_1TS,3))

# time step 1
meshViewer.ViewTime = times[1]
meshViewer.Representations[0].UpdatePipeline()
Render()
cr.UpdatePipeline()
ds0 = servermanager.Fetch(cr)
ts1 = numpy_support.vtk_to_numpy(ds0.GetPointData().GetAbstractArray(0))[::step]
MyAssert(not np.all(ts0==ts1))
MyAssert(ts1.shape == (NB_OF_LINES_FOR_1TS,3))

# time step 2
meshViewer.ViewTime = times[2]
meshViewer.Representations[0].UpdatePipeline()
Render()
cr.UpdatePipeline()
ds0 = servermanager.Fetch(cr)
ts2 = numpy_support.vtk_to_numpy(ds0.GetPointData().GetAbstractArray(0))[::step]
MyAssert(ts2.shape == (NB_OF_LINES_FOR_1TS-1,3)) # ???? for TS == 2.0 Plot AS is missing !?
