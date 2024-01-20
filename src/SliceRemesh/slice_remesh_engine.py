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
 
import medcoupling as mc
import vtk2medcoupling
import vtk
import tempfile
import os

def Remesh(in_mesh_path, out_mesh_path, constant_size):
    import salome
    salome.salome_init()
    import SMESH, SALOMEDS
    from salome.smesh import smeshBuilder
    smesh = smeshBuilder.New()

    #--- 1. Création de l'hypothèse ---
    hypo = smesh.CreateAdaptationHypothesis()
    #-- Ajout des options
    hypo.setSizeMapType('Constant')
    hypo.setConstantSize(constant_size)
    # Options avancées
    hypo.setOptionValue("adaptation", "surface")

    #--- 2. Création de l'objet pour l'adaptation ---
    objet_adapt = smesh.Adaptation('MG_Adapt')

    #--- Fichiers IN/OUT ---#
    objet_adapt.setMEDFileIn(in_mesh_path)
    objet_adapt.setMEDFileOut(out_mesh_path)

    #--- Association de l'hypothese à l'objet de l'adaptation
    objet_adapt.AddHypothesis(hypo)

    #print (hypo.getCommandToRun())

    #-- 3. Calcul - True pour publication dans SMESH
    err = objet_adapt.Compute(True)

    if err != 0:
        raise RuntimeError( "Error during computation of mesh !" )

def ReavaluateField(curr_field, new_mesh, eps = 1e-12):
    """
        Compute a curr field mesh on a new mesh
    """
    # Recompute the field on the extracted slice
    curr_field_mesh = curr_field.getMesh()
    loc = curr_field.getDiscretization().getLocalizationOfDiscValues(new_mesh)
    res = curr_field.getDiscretization().getValueOnMulti(curr_field.getArray(), curr_field_mesh, loc)
    ret = mc.MEDCouplingFieldDouble.New(curr_field.getDiscretization().getEnum ())
    ret.setArray(res) ; ret.setMesh(new_mesh) ; ret.copyAllTinyAttrFrom(curr_field)
    return ret

def FilterOnlyTopLevelCells( ugin2d ):
    import vtk
    from vtk.util import numpy_support
    ecbt = vtk.vtkExtractCellsByType()
    ecbt.SetInputData( ugin2d )
    allTypes = numpy_support.vtk_to_numpy( ugin2d.GetCellTypesArray() )

    sel0 = [elt.split("_")[1] for elt in dir(vtk) if "VTK_"==elt[:4] and len(elt.split("_"))==2]
    sel1 = ["{}{}".format(elt[0],elt[1:].lower()) for elt in sel0]
    sel2 = [elt for elt in sel1 if "vtk{}".format(elt) in dir(vtk)]
    sel3 = []
    for elt in sel2:
        try:
            inst = eval("vtk.vtk{}()".format(elt))
            inst.clsname = elt
            sel3.append( inst )
        except:
            pass
    sel4 = [elt for elt in sel3 if isinstance(elt,vtk.vtkCell)]
    sel5 = {elt.GetCellType() : elt for elt in sel4}
    if len(sel4) != len(sel5):
        raise RuntimeError("Presence of multiple sub vtkCell instance having same CellType !")
    dimToSelect = max([sel5[elt].GetCellDimension() for elt in allTypes])
    ecbt.RemoveAllCellTypes()
    for elt in [sel5[elt].GetCellType() for elt in set(allTypes)]:
        ecbt.AddCellType(elt)
    ecbt.Update()
    ugin = ecbt.GetOutputDataObject(0)
    return ugin

def ConvertUGToMCFields(ugin):
    mc_mesh_in = vtk2medcoupling.mesh_convertor_mem(ugin)
    # on ne prend que les champs float64
    dblsArr = [ugin.GetPointData().GetArray(i) for i in range(ugin.GetPointData().GetNumberOfArrays()) if isinstance( ugin.GetPointData().GetArray(i), vtk.vtkDoubleArray )]
    from vtk.util import numpy_support
    mc_mesh_in.unPolyze()
    o2n = mc_mesh_in.rearrange2ConsecutiveCellTypes()
    # on ne prend que le 1er champ
    ret = []
    for dblArr in dblsArr:
        arr = mc.DataArrayDouble(numpy_support.vtk_to_numpy(dblArr))
        arr = arr.renumber(o2n)
        f = mc.MEDCouplingFieldDouble(mc.ON_NODES)
        f.setArray( arr )
        f.setMesh( mc_mesh_in )
        f.setName( dblArr.GetName() )
        ret.append( f )
    return ret

def ConvertMCFieldToVTKUG( mcfields ):
    """
    :param mcfield: list of fields to be converted to unstructured grid
    :type mcfield: list of mc.MEDCouplingFieldDouble
    """
    if len(mcfields) < 1:
        raise RuntimeError("Input list must be of size >= 1 !")
    def ConvertMCFieldToVTKUGInternal( mcfield ):
        with tempfile.TemporaryDirectory() as d:
            tmpFileName = os.path.join(d,"resu.vtu")
            effectiveFileName = mcfield.writeVTK(tmpFileName)
            if effectiveFileName != tmpFileName:
                raise RuntimeError("input MEDCoupling Field seems to lie on invalid of non unstructuredgrid")
            ugr = vtk.vtkXMLUnstructuredGridReader()
            ugr.SetFileName(tmpFileName)
            ugr.Update()
            resu = ugr.GetOutputDataObject(0)
        return resu
    dsRet = ConvertMCFieldToVTKUGInternal( mcfields[0] )
    for mcfield in mcfields[1:]:
        ds = ConvertMCFieldToVTKUGInternal( mcfield )
        dsRet.GetPointData().AddArray( ds.GetPointData().GetArray(0) )
    return dsRet

def EngineOfRemesh(ugin2d,ugin3d,constantSize,mergeNodesTol):
    """
    :param ugin2d: 2D dataset that will be remeshed. Arrays included in ugin2d are simply ignored only cells/points are considered.
    :type ugin2d: vtk.vtkUnstructuredGrid
    :param ugin3d: 3D dataset used to resample values in the remeshed 2D mesh ugin2d.
    :param constantSize: Float indicating the size of cells to be generated
    :param mergeNodesTol: Float indicating the tolerance below which nodes in ugin2d will be merged.
    """
    uginclean2d = FilterOnlyTopLevelCells( ugin2d )
    mesh2d_mc = vtk2medcoupling.mesh_convertor_mem(uginclean2d)
    mesh2d_mc.unPolyze()
    mesh2d_mc.mergeNodes(mergeNodesTol)
    uginclean3d = FilterOnlyTopLevelCells( ugin3d )
    with tempfile.TemporaryDirectory() as d:
        tempFileNameIn = os.path.join(d,"inpp.med")
        tempFileNameOut = os.path.join(d,"outpp.med")
        mesh2d_mc.write( tempFileNameIn )
        Remesh(tempFileNameIn,tempFileNameOut,constantSize)
        tempFileCvt = os.path.join(d,"mesh.vtu")
        mc.MEDFileMesh.New(tempFileNameOut)[0].writeVTK(tempFileCvt)
        rd1 = vtk.vtkXMLUnstructuredGridReader()
        rd1.SetFileName(tempFileCvt)
        rd1.Update()
        #
        rds = vtk.vtkResampleWithDataSet()
        rds.SetInputData(rd1.GetOutputDataObject(0))
        rds.SetSourceData(uginclean3d)
        rds.Update()
        return rds.GetOutputDataObject(0)

def EngineOfRemeshSubProcess(ugin2d,ugin3d,constantSize,mergeNodesTol):
    with tempfile.TemporaryDirectory() as d:
        import subprocess as sp
        ugw = vtk.vtkXMLUnstructuredGridWriter()
        ugw.SetInputData(ugin2d)
        par1 = os.path.join(d,"twod.vtu")
        ugw.SetFileName(par1)
        ugw.Update()
        #
        ugw = vtk.vtkXMLUnstructuredGridWriter()
        ugw.SetInputData(ugin3d)
        par2 = os.path.join(d,"three.vtu")
        ugw.SetFileName(par2)
        ugw.Update()
        #
        outputFileName = os.path.join(d,"return.vtu")
        cmd = ["python3",__file__,"-d",par1,"-t",par2,"-o",outputFileName,"-s",str(constantSize),"-z",str(mergeNodesTol)]
        print("Cmd launched : {}".format(" ".join(cmd)))
        p = sp.Popen(cmd)
        p.communicate()
        #
        rd1 = vtk.vtkXMLUnstructuredGridReader()
        rd1.SetFileName(outputFileName)
        rd1.Update()
        return rd1.GetOutputDataObject(0)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--twod', dest="twod", required=True, help="2D VTU dataset to remesh")
    parser.add_argument('-t','--threed', dest="threed", required=True, help="3D VTU dataset")
    parser.add_argument('-o','--output', dest="output", required=True, help="3D VTU dataset")
    parser.add_argument('-s','--constant-size', dest = "constantSize", type=float, default=0.04, help="Float indicating the size of cells to be generated")
    parser.add_argument('-z','--mergenodes-tolerance', dest = "mergeNodesTolerance", type=float, default=1e-7, help="Float indicating the tolerance below which nodes in ugin2d will be merged")
    args = parser.parse_args()
    rd1 = vtk.vtkXMLUnstructuredGridReader()
    rd1.SetFileName(args.twod)
    rd1.Update()
    rd2 = vtk.vtkXMLUnstructuredGridReader()
    rd2.SetFileName(args.threed)
    rd2.Update()
    #
    ret = EngineOfRemesh( rd1.GetOutputDataObject(0), rd2.GetOutputDataObject(0),  args.constantSize, args.mergeNodesTolerance)
    #
    ugw = vtk.vtkXMLUnstructuredGridWriter()
    ugw.SetInputData(ret)
    ugw.SetFileName(args.output)
    ugw.Update()
    pass
