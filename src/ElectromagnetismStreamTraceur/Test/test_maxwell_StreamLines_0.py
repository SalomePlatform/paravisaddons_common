from paraview.simple import *
from medcoupling import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

def MyAssert(clue):
    if not clue:
        raise RuntimeError("Assertion failed !")

fname = "maxwell_streamline.med"

arr = DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10])
m = MEDCouplingCMesh()
m.setCoords(arr,arr,arr)
m = m.buildUnstructured()

m.changeSpaceDimension(3,0.)
m.setName("mesh")
f = MEDCouplingFieldDouble(ON_CELLS)
f.setMesh(m)
f.setName("field")
arrf = DataArrayDouble(10*10*10,3)
arrf[:,0] = 1 ; arrf[:,1] = 0 ; arrf[:,2] = 0
f.setArray( arrf )
f.getArray().setInfoOnComponents(["X","Y","Z"])
f.checkConsistencyLight()
f.write(fname)
f2 = f.deepCopy()
arrf2 = DataArrayDouble(10*10*10,3)
arrf2[:,0] = 0 ; arrf2[:,1] = 1 ; arrf2[:,2] = 0
f2.setArray( arrf2 )
f2.setName("field2")
WriteFieldUsingAlreadyWrittenMesh(fname,f2)

testmed = MEDReader(FileNames=fname)
testmed.FieldsStatus = ['TS0/mesh/ComSup0/field@@][@@P0', 'TS0/mesh/ComSup0/field2@@][@@P0', 'TS0/mesh/ComSup0/mesh@@][@@P0']
testmed.TimesFlagsStatus = ['0000']
streamTraceur1 = LigneDeChamp(Input=testmed,SeedType='Point Cloud')
streamTraceur1.SeedType.Radius = 1
streamTraceur1.SeedType.Center = [ 7.23,7.26,3.42 ]
streamTraceur1.Vectors = ['CELLS', "field"]
streamTraceur1.UpdatePipeline()
ds0 = servermanager.Fetch(streamTraceur1)
MyAssert(ds0.GetNumberOfCells()==200)
