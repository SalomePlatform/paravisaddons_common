# Copyright (C) 2021-2024  CEA, EDF
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

CMAKE_POLICY(SET CMP0071 NEW)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Common CMake macros
# ===================
SET(TMP_CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})
unset(CMAKE_MODULE_PATH)
SET(CONFIGURATION_ROOT_DIR $ENV{CONFIGURATION_ROOT_DIR} CACHE PATH "Path to the Salome CMake configuration files")
IF(EXISTS ${CONFIGURATION_ROOT_DIR})
  LIST(APPEND CMAKE_MODULE_PATH "${CONFIGURATION_ROOT_DIR}/cmake")
  INCLUDE(SalomeMacros)
ELSE()
  MESSAGE(FATAL_ERROR "We absolutely need the Salome CMake configuration files, please define CONFIGURATION_ROOT_DIR !")
ENDIF()

SET(MEDCOUPLING_ROOT_DIR $ENV{MEDCOUPLING_ROOT_DIR} CACHE PATH "Path to the MEDCoupling tool")
IF(EXISTS ${MEDCOUPLING_ROOT_DIR})
  LIST(APPEND CMAKE_MODULE_PATH "${MEDCOUPLING_ROOT_DIR}/cmake_files")
ENDIF()
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_ROOT}/Modules")
LIST(APPEND CMAKE_MODULE_PATH ${TMP_CMAKE_MODULE_PATH})

INCLUDE(SalomeSetupPlatform)
SET(BUILD_SHARED_LIBS TRUE)

FIND_PACKAGE(SalomeHDF5 REQUIRED)
FIND_PACKAGE(SalomeMEDCoupling REQUIRED)

SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH NOCHECK ${CMAKE_INSTALL_PREFIX}/${SALOME_INSTALL_BINS}
                                                 ${CMAKE_INSTALL_PREFIX}/${SALOME_INSTALL_PYTHON})
SALOME_ACCUMULATE_ENVIRONMENT(LD_LIBRARY_PATH NOCHECK ${CMAKE_INSTALL_PREFIX}/${SALOME_INSTALL_LIBS})
SALOME_ACCUMULATE_ENVIRONMENT(PV_PLUGIN_PATH NOCHECK ${CMAKE_INSTALL_PREFIX}/lib/paraview)

paraview_add_plugin(ZJFilterPlugin
  VERSION "1.0"
  MODULES ZJFilterModule
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/ZJFilterModule/vtk.module"
  SERVER_MANAGER_XML filters.xml
)

install(TARGETS ZJFilterPlugin
  RUNTIME DESTINATION lib/paraview
  LIBRARY DESTINATION lib/paraview
  ARCHIVE DESTINATION lib/paraview
)
