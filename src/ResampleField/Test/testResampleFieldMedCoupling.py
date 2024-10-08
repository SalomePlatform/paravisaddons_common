from paraview.simple import *

import vtk
import medcoupling as mc
from vtk.util import numpy_support
import numpy as np

from pathlib import Path
import tempfile

arr = mc.DataArrayDouble(3) ; arr.iota()
mSrc = mc.MEDCouplingCMesh() ; mSrc.setCoords(arr,arr) ; mSrc = mSrc.buildUnstructured() ; mSrc.changeSpaceDimension(3,0.)
fCellSrc = mc.MEDCouplingFieldDouble(mc.ON_CELLS)
fCellSrc.setMesh( mSrc )
fCellSrc.setName("CellField")
fCellSrc.setArray( mc.DataArrayDouble( [0,1,2,3] ) )
mTrg = mSrc.deepCopy()
mTrg.renumberNodes([5,3,2,1,0,8,7,6,4],mTrg.getNumberOfNodes())
# srcCell 0 goes to position 3 in targetMesh
# srcCell 1 goes to position 2 in targetMesh
# srcCell 2 goes to position 0 in targetMesh
# srcCell 3 goes to position 1 in targetMesh
mTrg.renumberCells([3,2,0,1])

with tempfile.TemporaryDirectory() as di:
    srcFname = Path(di) / "src.vtu"
    mc.MEDCouplingFieldDouble.WriteVTK(f"{srcFname}", [fCellSrc] )
    srcvtu = XMLUnstructuredGridReader(FileName=[f"{srcFname}"])
    srcvtu.CellArrayStatus = [fCellSrc.getName()]
    trgFname = Path(di) / "trg.vtu"
    mTrg.writeVTK( f"{trgFname}" )
    trgvtu = XMLUnstructuredGridReader(FileName=[f"{trgFname}"])
    trgvtu.CellArrayStatus = []
    resampleField1 = ResampleField( SourceData=srcvtu, TargetMesh=trgvtu )
    resampleField1.UpdatePipeline()
    res = servermanager.Fetch( resampleField1 )
    resValueNp = numpy_support.vtk_to_numpy( res.GetAttributes(vtk.VTK_CELL_MODE).GetArray("CellField") )
    ref = np.array([2., 3., 1., 0.], dtype=np.float32)
    assert( np.all( np.isclose( resValueNp, ref, rtol = 0., atol = 1e-5 ) ) )
