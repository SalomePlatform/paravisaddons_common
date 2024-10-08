"""A Python Plugin for ParaView which resample MED data to have a the same order"""

from paraview.util.vtkAlgorithm import *

import vtk
# TODO: n'ayant pas accès à ces 2 librairies je ne peux être sur que cela soit suffisant
# il faudra vous assurez que pvpython connaissent ces 2 librairies
import medcoupling as mc
from vtk2medcoupling import mesh_convertor_mem
import numpy as np

#------------------------------------------------------------------------------
def resampleMedCoupling(srcVtk, trgVtk, eps):
    meshTrgMed = mesh_convertor_mem( trgVtk )
    meshSrcMed = mesh_convertor_mem( srcVtk )
    dataffvtk = srcVtk.GetAttributes(vtk.VTK_CELL_MODE)
    srcFieldsMed = [ vtkArray_to_mc_field( dataffvtk.GetArray(i),meshSrcMed  ) for i in range( dataffvtk.GetNumberOfArrays() ) ]
    fieldsMC = remapSimple( srcFieldsMed, meshTrgMed, eps)
    from pathlib import Path
    import tempfile
    with tempfile.TemporaryDirectory() as di:
        fname = Path(di) / "temp.vtu"
        mc.MEDCouplingFieldDouble.WriteVTK( f"{fname}", fieldsMC )
        rd = vtk.vtkXMLUnstructuredGridReader()
        rd.SetFileName(fname)
        rd.Update()
        return rd.GetOutputDataObject(0)

#------------------------------------------------------------------------------
def remapSimple(srcFields, targetMesh, eps):
    srcMesh = srcFields[0].getMesh()
    res = targetMesh.checkGeoEquivalWith(srcMesh,12,eps)
    ret = []
    for srcField in srcFields:
        ff = srcField.deepCopy()
        ff.renumberCells(res[0])
        ff.setMesh(targetMesh)
        ret.append( ff )
    return ret
    
#------------------------------------------------------------------------------
def NumpyToMC64(nparr):
    return mc.DataArrayDouble(nparr)

#------------------------------------------------------------------------------
def NumpyToMC32(nparr):
    return mc.DataArrayFloat(nparr).convertToDblArr()

#------------------------------------------------------------------------------
def NumpyToMCInt32(nparr):
    return mc.DataArrayInt32(nparr).convertToDblArr()

Dico = {np.float32 : NumpyToMC32, np.float64: NumpyToMC64, np.int32 : NumpyToMCInt32}

#------------------------------------------------------------------------------
def vtkArray_to_mc_field( vtkarr, mcMesh ):
    def CompoNameSafe(compName):
        if compName is None:
            return ""
        else:
            return compName
    from vtk.util import numpy_support
    f = mc.MEDCouplingFieldDouble(mc.ON_CELLS)
    f.setName( vtkarr.GetName() )
    nparr = numpy_support.vtk_to_numpy( vtkarr )
    f.setArray( Dico[nparr.dtype.type](nparr) )
    f.setMesh(mcMesh)
    f.getArray().setInfoOnComponents( [CompoNameSafe( vtkarr.GetComponentName(i) ) for i in range( vtkarr.GetNumberOfComponents() ) ] )
    return f

#------------------------------------------------------------------------------
@smproxy.filter(name="Resample Field")
#@smproperty.input(name="InputUnstructuredGrid", port_index=0)
#@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
@smproperty.xml("""
<InputProperty command="SetInputConnection" name="SourceData" port_index="0">
  <ProxyGroupDomain name="groups">
    <Group name="sources" />
    <Group name="filters" />
  </ProxyGroupDomain>
  <DataTypeDomain name="input_type">
    <DataType value="vtkUnstructuredGrid" />
  </DataTypeDomain>
</InputProperty>
<InputProperty command="SetInputConnection" name="TargetMesh" port_index="1">
  <ProxyGroupDomain name="groups">
    <Group name="sources" />
    <Group name="filters" />
  </ProxyGroupDomain>
  <DataTypeDomain name="input_type">
    <DataType value="vtkUnstructuredGrid" />
  </DataTypeDomain>
</InputProperty>
""")
class ResampleFieldUsingMEDCoupling(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=2, nOutputPorts=1, outputType="vtkUnstructuredGrid")

    def FillInputPortInformation(self, port, info):
        if port == 0:
            info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        srcData = self.GetInputData(inInfoVec, 0, 0)
        trgMesh = self.GetInputData(inInfoVec, 1, 0)
        outData = self.GetOutputData(outInfoVec, 0)

        #from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
        #input = vtkUnstructuredGrid.GetData(inInfoVec[0], 0)
        #output = vtkUnstructuredGrid.GetData(outInfoVec, 0)
        
        output = resampleMedCoupling(srcData,trgMesh,1e-5)
        outData.ShallowCopy(output)

        return 1
