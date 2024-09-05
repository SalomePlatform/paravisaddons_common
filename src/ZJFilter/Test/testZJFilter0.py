from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
from vtk.util import numpy_support

def MyAssert(clue):
    if not clue:
        raise RuntimeError("Assertion failed !")

fileName="testZJFilter.med"
testmed = MEDReader(FileNames=fileName)
testmed.FieldsStatus = ['TS0/CONTACTOR_BTBrin_BT1_USI/ComSup0/pertes_Ohm@@][@@P0']
testmed.TimesFlagsStatus = ['0000']

# create a new 'ZJ Filter'
zJFilter1 = ZJFilter(registrationName='ZJFilter1', Input=testmed)
zJFilter1.UpdatePipeline()
cdi = zJFilter1.GetCellDataInformation()
MyAssert(cdi.GetNumberOfArrays() == 5)
print("OK")