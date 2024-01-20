#  -*- coding: utf-8 -*-
# Copyright (C) 2022-2024  CEA, EDF
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

from paraview.util.vtkAlgorithm import *
import importlib.util
import sys
import os
from paraview import print_info, print_error
from vtkmodules.vtkCommonDataModel import vtkCompositeDataSet
import vtk
import tempfile
import os

@smproxy.filter(name="Slice Remesh")
@smproperty.xml("""
<InputProperty command="SetInputConnection" name="Input3d" port_index="0">
  <ProxyGroupDomain name="groups">
    <Group name="sources" />
    <Group name="filters" />
  </ProxyGroupDomain>
  <DataTypeDomain name="input_type">
    <DataType value="vtkDataSet" />
    <DataType value="vtkCompositeDataSet" />
  </DataTypeDomain>
</InputProperty>
<InputProperty command="SetInputConnection" name="InputSliceToRemesh2d" port_index="1">
  <ProxyGroupDomain name="groups">
    <Group name="sources" />
    <Group name="filters" />
  </ProxyGroupDomain>
  <DataTypeDomain name="input_type">
    <DataType value="vtkDataSet" />
    <DataType value="vtkCompositeDataSet" />
  </DataTypeDomain>
</InputProperty>
""")
class SliceRemesh(VTKPythonAlgorithmBase):
  """
  This filter allows to hotplug a python module at runtime to process an input.
  The python module must respect these conditions :
    - depend on package that will be available in ParaView Python environement ONLY
    - contain a function `remesh()` that returns a vtkUnstructuredGrid instance and accept 2 arguments
      * First argument is the input dataset
      * Second argument is a dictionary of parameters to pass to the mesher. Currently exposed key are : "target_size".

  If the environment variable 'SLICERMESH_DEFAULT_MESHER' exists and points to a valid file then this file will be used
  as default module.
  """

  MODULE_NAME = "PythonMesher.mesher"

  def __init__(self):
    super().__init__(nInputPorts=2, nOutputPorts=1, outputType="vtkDataObject")
    self._mesherParameters = {
      "target_size": 1.0,
      "tolerance_merge_nodes":1.0
    }

  def FillInputPortInformation(self, port, info):
    if port == 0 or port == 1:
      info.Remove(self.INPUT_REQUIRED_DATA_TYPE())
      info.Append(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet")
      info.Append(self.INPUT_REQUIRED_DATA_TYPE(), "vtkCompositeDataSet")
      return 1
    else:
      return 0

  def RequestDataObject(self, request, inInfo, outInfo):
    inData = self.GetInputData(inInfo, 0, 0)
    outData = self.GetOutputData(outInfo, 0)
    assert inData is not None

    outData = vtk.vtkUnstructuredGrid()
    outInfo.GetInformationObject(0).Set(outData.DATA_OBJECT(), outData)
    return super().RequestDataObject(request, inInfo, outInfo)

  @smproperty.doublevector(name="TargetSize", label="Target Cell Size", number_of_elements=1, default_values=0.04)
  def SetMeshSize(self, size):
    """
    Set the target size of cells for the remesher.
    """
    if self._mesherParameters["target_size"] != size:
      self._mesherParameters["target_size"] = size
      self.Modified()

  @smproperty.doublevector(name="ToleranceMergeNodes", label="Tolerance Merge Nodes", number_of_elements=1, default_values=1e-7)
  def SetToleranceMergeNodes(self, size):
    """
    Set the tolerance use to merge nodes.
    """
    if self._mesherParameters["tolerance_merge_nodes"] != size:
      self._mesherParameters["tolerance_merge_nodes"] = size
      self.Modified()

  def RequestData(self, request, inInfo, outInfo):
    in3dData = self.GetInputData(inInfo, 0, 0)
    inSliceToRemesh2dData = self.GetInputData(inInfo, 1, 0)
    outData = self.GetOutputData(outInfo, 0)

    if (inSliceToRemesh2dData.IsA("vtkCompositeDataSet")):
      inSliceToRemesh2dData = inSliceToRemesh2dData.GetBlock(0)

    if (in3dData.IsA("vtkCompositeDataSet")):
      in3dData = in3dData.GetBlock(0)

    if inSliceToRemesh2dData.IsA("vtkPolyData"):
      # vtk.vtkAppendFilter
      ret = vtk.vtkUnstructuredGrid()
      ret.DeepCopy(inSliceToRemesh2dData)
      inSliceToRemesh2dData = ret ; del ret

    import slice_remesh_engine
    ret = slice_remesh_engine.EngineOfRemeshSubProcess(inSliceToRemesh2dData,in3dData,self._mesherParameters["target_size"],self._mesherParameters["tolerance_merge_nodes"])
    outData.ShallowCopy(ret)

    return 1
