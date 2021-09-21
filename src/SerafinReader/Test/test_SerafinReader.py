from paraview.simple import *
from vtk.util import numpy_support
import numpy as np

paraview.simple._DisableFirstRenderCameraReset()

# Normal single precision file
reader = SerafinReader(FileName='geo_sp.slf')
reader.UpdatePipeline()
reader_ds = servermanager.Fetch(reader)
# Checking that fields THING X THING Y are merged into THING *
list_var = [reader_ds.GetPointData().GetArrayName(i).strip() for i in range(reader_ds.GetPointData().GetNumberOfArrays())]
assert("THING *" in list_var)
assert("THING X" not in list_var)
assert("THING OTHER" in list_var)
# Checking that fields VELOCITY U, VELOCITY V are merged into VELOCITY *
assert("VELOCITY *" in list_var)
assert("VELOCITY U" not in list_var)
assert("VELOCITY STUFF" in list_var)
# Checking that VELOCITY has 3 components shoulb be on 0
assert(reader_ds.GetPointData().GetArray(0).GetNumberOfComponents() == 3)

# Simple read of Litte Endian file (Would crash if issue in the plugin)
reader = SerafinReader(FileName='geo_sp_le.slf')
# Simple read of Litte Endian file (Would crash if issue in the plugin)
reader = SerafinReader(FileName='geo_dp_le.slf')

# Normal double precision file
reader = SerafinReader(FileName='geo_dp.slf')
reader.UpdatePipeline()
reader_ds = servermanager.Fetch(reader)
# Checking that values are properly read
data = numpy_support.vtk_to_numpy(reader_ds.GetPointData().GetArray(0))
assert(np.all(data <= 1e-8))

# 3d mesh
reader = SerafinReader(FileName='geo_3d_dp.slf')
reader.UpdatePipeline()
reader_ds = servermanager.Fetch(reader)

# Reading file with time steps (10)
reader = SerafinReader(FileName='result.slf')
reader.UpdatePipeline()
reader_ds = servermanager.Fetch(reader)

# Using animationScene to get tim information
anim = GetAnimationScene()
tk = anim.TimeKeeper
times = tk.TimestepValues
ref_times = np.array([i*1.0 for i in range(10)])
diff = np.abs(ref_times - times)
assert(np.all(diff <= 1e-6))
