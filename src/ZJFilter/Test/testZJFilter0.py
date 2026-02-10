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
ref_set = {'FamilyIdCell', 'NumIdCell', 'BT1Cond01', 'pertes_Ohm'}
set_to_test = set( [cdi.GetArray(i).GetName() for i in range( cdi.GetNumberOfArrays()  ) ] )
print( set_to_test )
MyAssert(ref_set <= set_to_test)
print("OK")