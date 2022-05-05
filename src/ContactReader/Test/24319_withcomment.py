# Copyright (C) 2021-2022  CEA/DEN, EDF R&D
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

__doc__ = """ non regression test of EDF24319
resultante_cdth.rco -> 24319_no_comment.rco
resultante_cdth.resu -> 24319_with_comments.rco
"""

from paraview.simple import *

def MyAssert(b):
  if not b:
    raise RuntimeError("Oooops : Assertion failed !")


cr_without_comment = ContactReader(registrationName='WithoutComment', FileName='24319_no_comment.rco')
cr_with_comments = ContactReader(registrationName='WithComments', FileName='24319_with_comments.rco')
cr_without_comment.UpdatePipeline()
cr_with_comments.UpdatePipeline()

polydata_without_comment = servermanager.Fetch(cr_without_comment)
polydata_with_comments = servermanager.Fetch(cr_with_comments)
MyAssert(polydata_without_comment.GetNumberOfCells() == polydata_with_comments.GetNumberOfCells())
MyAssert(polydata_without_comment.GetNumberOfPoints() == polydata_with_comments.GetNumberOfPoints())
