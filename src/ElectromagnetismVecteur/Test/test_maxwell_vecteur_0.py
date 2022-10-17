from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
from vtk.util import numpy_support
import numpy as np
import medcoupling as mc

def MyAssert(clue):
    if not clue:
        raise RuntimeError("Assertion failed !")

fileName = "MaxwellVecteur.med"

# on genere un use case

arr = mc.DataArrayDouble([0,1,2])
m = mc.MEDCouplingCMesh()
m.setCoords(arr,arr)
m = m.buildUnstructured()
m.changeSpaceDimension(3,0.)
m.setName("mesh")
f = mc.MEDCouplingFieldDouble(mc.ON_CELLS)
f.setMesh(m)
f.setName("field")
f.setArray( mc.DataArrayDouble([(1,1,0), (2,-2,0), (-3,3,0), (-4,-4,0)]) )
f.getArray().setInfoOnComponents(["X","Y","Z"])
f.checkConsistencyLight()
f.write(fileName)
f2 = f.deepCopy()
f2.setArray( mc.DataArrayDouble([(1,0,0), (0,-2,0), (-3,0,0), (0,4,0)]) )
f2.setName("field2")
mc.WriteFieldUsingAlreadyWrittenMesh(fileName,f2)

# fin de la generation du fichier MED on attaque le test a proprement parler

testmed = MEDReader(FileNames=fileName)
testmed.FieldsStatus = ['TS0/mesh/ComSup0/field@@][@@P0', 'TS0/mesh/ComSup0/field2@@][@@P0', 'TS0/mesh/ComSup0/mesh@@][@@P0']
testmed.TimesFlagsStatus = ['0000']

vecteur1 = Vecteur(Input=testmed)
fieldName = "field"
colorFieldName = vecteur1.GetProperty("ColorArrayName")
vecteur1.OrientationArray = ['POINTS', fieldName]
vecteur1.ScaleFactor = 0.1
vecteur1.UpdatePipeline()

ds0 = servermanager.Fetch(vecteur1)
MyAssert(ds0.GetNumberOfBlocks()==1)
ds0 = ds0.GetBlock(0)

arr = numpy_support.vtk_to_numpy(ds0.GetPointData().GetArray(fieldName))
arr2 = numpy_support.vtk_to_numpy(ds0.GetPointData().GetArray(fieldName))
MyAssert( np.all(np.isclose(arr,arr2,0,1e-12)) ) # on check que la couleur est aligné sur le field selectionné
MyAssert( np.all(np.isclose(arr[-1],np.array([-4.,-4.,0.]),0,1e-12)) )
MyAssert( arr.shape==(20,3) ) # 20 == 4 * 5 . Une fleche c est 5 points (3 pour la tete et 2 pour la tige). 4 fleches car 4 cells

fieldName = "field2"
vecteur1.OrientationArray = ['POINTS', fieldName]
vecteur1.UpdatePipeline()
ds0 = servermanager.Fetch(vecteur1).GetBlock(0)
arr = numpy_support.vtk_to_numpy(ds0.GetPointData().GetArray(fieldName))
arr2 = numpy_support.vtk_to_numpy(ds0.GetPointData().GetArray(fieldName))
MyAssert( np.all(np.isclose(arr,arr2,0,1e-12)) )# on check que la couleur est aligné sur le field selectionné
MyAssert( np.all(np.isclose(arr[-1],np.array([0.,4.,0.]),0,1e-12)) )
MyAssert( arr.shape==(20,3) )
# on check la gueule des fleches notamment le fait que comme demandé on n'a pas de scale array avec des fleches de meme tailles bien qu avec des magnitudes differentes
pts = numpy_support.vtk_to_numpy(ds0.GetPoints().GetData())
pts_ref = np.array([[0.45, 0.5, 0], [0.55, 0.5, 0], [0.52, 0.49, 0], [0.55, 0.5, 0], [0.52, 0.51, 0], [1.5, 0.55, 0.], [1.5, 0.45, 0.], [1.51, 0.48, 0.], [1.5, 0.45, 0.], [1.49, 0.48, 0.], [0.55, 1.5, 0.], [0.45, 1.5, 0.],[0.48, 1.49, 0.], [0.45, 1.5, 0.], [0.48, 1.51, 0.], [1.5, 1.45, 0.], [1.5, 1.55, 0.], [1.49, 1.52, 0.], [1.5, 1.55, 0.], [1.51, 1.52, 0.]], dtype = np.float32)
MyAssert( np.all(np.isclose(pts,pts_ref,0,1e-12)) )
vecteur1.ScaleFactor = 1.0
vecteur1.UpdatePipeline()
ds0 = servermanager.Fetch(vecteur1).GetBlock(0)
pts = numpy_support.vtk_to_numpy(ds0.GetPoints().GetData())
pts_ref_1 = np.array([(0, 0.5, 0), (1, 0.5, 0), (0.7, 0.4, 0), (1, 0.5, 0), (0.7, 0.6, 0), (1.5, 1, 0), (1.5, 0, 0), (1.6, 0.3, 0), (1.5, 0, 0), (1.4, 0.3, 0), (1, 1.5, 0), (0, 1.5, 0), (0.3, 1.4, 0), (0, 1.5, 0), (0.3, 1.6, 0), (1.5, 1, 0), (1.5, 2, 0), (1.4, 1.7, 0), (1.5, 2, 0), (1.6, 1.7, 0)],dtype=np.float32)
