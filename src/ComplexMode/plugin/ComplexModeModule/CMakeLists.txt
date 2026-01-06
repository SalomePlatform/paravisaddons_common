# Copyright (C) 2021-2026  CEA, EDF
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

set(classes
  vtkComplexMode
)
# Visual C++ 2017 bug: defining _USE_MATH_DEFINES inside .cpp/.h file is inoperant
IF(WIN32)
ADD_COMPILE_DEFINITIONS(_USE_MATH_DEFINES)
ENDIF(WIN32)

vtk_module_add_module(ComplexModeModule
  FORCE_STATIC
  CLASSES ${classes}
  )
