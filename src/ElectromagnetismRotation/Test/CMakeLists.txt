# Copyright (C) 2012-2024  CEA, EDF
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
# Author : Yoann Audouin (EDF)

SET(test_maxwell_rotation_0_src
  mesh_benjamin_8_sept_2020.med
  test_maxwell_rotation_0.py
  )

SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)

ADD_TEST(ElectromagnetismRotation python test_maxwell_rotation_0.py)
SET_TESTS_PROPERTIES(ElectromagnetismRotation PROPERTIES ENVIRONMENT "${tests_env}")

# Application tests

SET(TEST_INSTALL_DIRECTORY ${PARAVISADDONS_COMMON_INSTALL_TESTS}/ElectromagnetismRotation)
INSTALL(FILES ${test_maxwell_rotation_0_src} DESTINATION ${TEST_INSTALL_DIRECTORY})

INSTALL(FILES CTestTestfileInstall.cmake
        DESTINATION ${TEST_INSTALL_DIRECTORY}
        RENAME CTestTestfile.cmake)
