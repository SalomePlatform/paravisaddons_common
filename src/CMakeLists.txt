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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CONFIGURATION_ROOT_DIR $ENV{CONFIGURATION_ROOT_DIR} CACHE PATH "Path to the Salome CMake files")
SET(PARAVISADDONS_COMMON_INSTALL_TESTS tests CACHE PATH "Install path: PARAVISADDONS_COMMON tests")
IF(EXISTS ${CONFIGURATION_ROOT_DIR})
  LIST(APPEND CMAKE_MODULE_PATH "${CONFIGURATION_ROOT_DIR}/cmake")
  INCLUDE(SalomeMacros)
ELSE()
  MESSAGE(FATAL_ERROR "We absolutely need the Salome CMake configuration files, please define CONFIGURATION_ROOT_DIR !")
ENDIF()

INCLUDE(SalomeSetupPlatform)

FOREACH(_Qt5_COMPONENT_ ${Qt5_FIND_COMPONENTS} ${Qt5_OPTIONAL_COMPONENTS})
  SET(_Qt5_COMPONENT Qt5${_Qt5_COMPONENT_})
  LIST(FIND Qt5_OPTIONAL_COMPONENTS ${_Qt5_COMPONENT_} idx)
  IF(${idx} GREATER -1)
    SET(Salome${_Qt5_COMPONENT}_FIND_QUIETLY TRUE)
  ENDIF()
  FIND_PACKAGE(${_Qt5_COMPONENT})
  LIST(APPEND QT_INCLUDES    ${${_Qt5_COMPONENT}_INCLUDE_DIRS})
  LIST(APPEND QT_DEFINITIONS ${${_Qt5_COMPONENT}_DEFINITIONS})
  LIST(APPEND QT_LIBRARIES   ${${_Qt5_COMPONENT}_LIBRARIES})
ENDFOREACH()

ENABLE_TESTING()

ADD_SUBDIRECTORY(MoveZCote)
ADD_SUBDIRECTORY(ComplexMode)
ADD_SUBDIRECTORY(QuadraticToLinear)
ADD_SUBDIRECTORY(ExtractComponentsPlugin)
ADD_SUBDIRECTORY(ZJFilter)
ADD_SUBDIRECTORY(SerafinReader)
ADD_SUBDIRECTORY(SpatialPfl)
ADD_SUBDIRECTORY(SinusXReader)
ADD_SUBDIRECTORY(TemporalOnPoint)
ADD_SUBDIRECTORY(DepthVsTime)
ADD_SUBDIRECTORY(SphereAlongLines)
ADD_SUBDIRECTORY(RateOfFlowThroughSection)
ADD_SUBDIRECTORY(ProbePointOverTime)
ADD_SUBDIRECTORY(CustomFilters)
ADD_SUBDIRECTORY(CellDataContour)
ADD_SUBDIRECTORY(AutoConvertPropertiesPlugin)
ADD_SUBDIRECTORY(AppendAttributesOverTime)
ADD_SUBDIRECTORY(ExtractThreeD)
ADD_SUBDIRECTORY(ContactReader)
ADD_SUBDIRECTORY(TorseurCIH)
ADD_SUBDIRECTORY(RosetteCIH)
ADD_SUBDIRECTORY(XYChartRepresentationColumns)
ADD_SUBDIRECTORY(GlyphCIH) # EDF15785
ADD_SUBDIRECTORY(ElectromagnetismVecteur)
ADD_SUBDIRECTORY(ElectromagnetismStreamTraceur)
ADD_SUBDIRECTORY(ElectromagnetismFluxDisc)
ADD_SUBDIRECTORY(ElectromagnetismRotation)


# Tests
configure_file(CTestTestfileInstall.cmake.in "CTestTestfileST.cmake" @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfileST.cmake DESTINATION ${PARAVISADDONS_COMMON_INSTALL_TESTS} RENAME CTestTestfile.cmake)
