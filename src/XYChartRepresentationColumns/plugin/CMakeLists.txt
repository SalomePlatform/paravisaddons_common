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

cmake_minimum_required(VERSION 3.8)
project(XYChartRepresentationColumnsPlugin)
find_package(ParaView REQUIRED)

set(BUILD_SHARED_LIBS ON)

paraview_add_plugin(XYChartRepresentationColumns
  VERSION "1.0"
  MODULES XYChartRepresentationColumnsModule
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/XYChartRepresentationColumnsModule/vtk.module"
  SERVER_MANAGER_XML views.xml
  XML_DOCUMENTATION OFF
)

set_target_properties (XYChartRepresentationColumns PROPERTIES POSITION_INDEPENDENT_CODE ON)
install(TARGETS XYChartRepresentationColumns
  RUNTIME DESTINATION lib/paraview
  LIBRARY DESTINATION lib/paraview
  ARCHIVE DESTINATION lib/paraview
)
