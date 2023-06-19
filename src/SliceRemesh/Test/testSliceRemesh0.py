#!/usr/bin/env python3
#  -*- coding: utf-8 -*-
# Copyright (C) 2020-2023  CEA, EDF
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

import sys
import salome

salome.salome_init()

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
Cylinder_1 = geompy.MakeCylinderRH(100, 300)

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

Mesh_1 = smesh.Mesh(Cylinder_1,'Mesh_1')
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Simple_Parameters_1 = NETGEN_1D_2D_3D.Parameters(smeshBuilder.SIMPLE)
NETGEN_3D_Simple_Parameters_1.SetNumberOfSegments( 15 )
NETGEN_3D_Simple_Parameters_1.LengthFromEdges()
NETGEN_3D_Simple_Parameters_1.LengthFromFaces()
isDone = Mesh_1.Compute()

if not isDone:
  raise RuntimeError("Fail to compute mesh !")

## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_3D_Simple_Parameters_1, 'NETGEN 3D Simple Parameters_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


fname = "toto.vtu"
mm = Mesh_1.ExportMEDCoupling()
import medcoupling as mc
import vtk
f = mc.MEDCouplingFieldDouble(mc.ON_NODES)
f.setMesh(mm[0])
arr = (f.getMesh().getCoords()-[0,0,150]).magnitude()
f.setArray(arr)
f.setName("FieldNode")
import tempfile
import os
with tempfile.TemporaryDirectory() as d:
  fullFileName = os.path.join(d,fname)
  f.writeVTK(fullFileName)

  ur = vtk.vtkXMLUnstructuredGridReader()
  ur.SetFileName(fullFileName)
  ur.Update()
  ur.GetOutput()


  from paraview.simple import *
  producer = TrivialProducer()
  producer.GetClientSideObject().SetOutput(ur.GetOutput())
  producer.UpdatePipeline()

  slice1 = Slice(Input=producer)
  slice1.SliceType = 'Plane'
  slice1.HyperTreeGridSlicer = 'Plane'
  slice1.SliceOffsetValues = [0.0]
  slice1.SliceType.Origin = [0.0, 0., 150.0]
  slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, 150.0]
  slice1.SliceType.Normal = [0.8560156765136049, -0.3977749519024914, 0.33017003074465473]
  slice1.UpdatePipeline()

  vtp = servermanager.Fetch(slice1)
  ret = vtk.vtkUnstructuredGrid()
  ret.DeepCopy(vtp)

  sliceUnstructured = TrivialProducer()
  sliceUnstructured.GetClientSideObject().SetOutput(ret)
  sliceUnstructured.UpdatePipeline()

  # test avec unstructured en input
  sliceRemesh1 = SliceRemesh(Input3d=producer,InputSliceToRemesh2d=sliceUnstructured)
  sliceRemesh1.TargetCellSize = 0.04
  sliceRemesh1.UpdatePipeline()
  res0 = servermanager.Fetch(sliceRemesh1)
  if res0.GetNumberOfCells() < 1000:
    raise RuntimeError("Test fails #1 !")

  sliceRemesh1.TargetCellSize = 1.0
  sliceRemesh1.UpdatePipeline()
  res1 = servermanager.Fetch(sliceRemesh1)
  if res1.GetNumberOfCells() < 1000 and res1.GetNumberOfCells() > 10:
    raise RuntimeError("Test fails #2 !")

  if res0.GetNumberOfCells() < res1.GetNumberOfCells() and res1.GetNumberOfCells() > 10:
    raise RuntimeError("Test fails #3 !")

  # test avec polydata en input pour le slice
  sliceRemesh2 = SliceRemesh(Input3d=producer,InputSliceToRemesh2d=slice1)
  sliceRemesh2.TargetCellSize = 20
  sliceRemesh2.UpdatePipeline()
  res2 = servermanager.Fetch(sliceRemesh2)

  if res1.GetNumberOfCells() < res2.GetNumberOfCells() and res2.GetNumberOfCells() > 10:
    raise RuntimeError("Test fails #4 !")
