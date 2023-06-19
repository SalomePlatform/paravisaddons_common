# Copyright (C) 2021-2023  CEA, EDF
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

set(interfaces)
set(sources)

cmake_policy(SET CMP0071 OLD) # bug in ParaViewPlugin.cmake?

if(PARAVIEW_USE_QT)

  set(sources
    pqElectroRotationAbstractFieldsWidget.cxx
    pqElectroRotationGroupWidget.cxx)

  paraview_plugin_add_property_widget(
    KIND WIDGET
    TYPE "ElectroRotationGroupWidgetType"
    CLASS_NAME pqElectroRotationGroupWidget
    INTERFACES property_interfaces
    SOURCES property_sources)
  list(APPEND interfaces
    ${property_interfaces})
  list(APPEND sources
    ${property_sources})

endif(PARAVIEW_USE_QT)

paraview_add_plugin(ElectromagnetismRotation
  VERSION "5.0"
  UI_INTERFACES ${interfaces}
  SOURCES ${sources}
  UI_RESOURCES Resources/pqElectromagnetismRotation.qrc
  SERVER_MANAGER_XML Resources/ElectromagnetismRotation.xml
  MODULES ElectromagnetismRotationIO
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../ElectromagnetismRotationIO/vtk.module"
  )

install(TARGETS ElectromagnetismRotation
  RUNTIME DESTINATION lib/paraview
  LIBRARY DESTINATION lib/paraview
  ARCHIVE DESTINATION lib/paraview
)
