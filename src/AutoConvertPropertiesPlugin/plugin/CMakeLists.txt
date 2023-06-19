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

# Create an auto-start plugin. Auto start plugins provide callbacks that get
# called when the plugin is loaded and when the application shutsdown.
paraview_plugin_add_auto_start(
  CLASS_NAME pqAutoConvertPropertiesStarter
  STARTUP onStartup
  SHUTDOWN onShutdown
  INTERFACES autostart_interface
  SOURCES autostart_sources
)

# create a plugin for this starter
paraview_add_plugin(AutoConvertProperties
  VERSION "1.0"
  UI_INTERFACES ${autostart_interface}
  SOURCES ${autostart_sources} pqAutoConvertPropertiesStarter.cxx
)

target_include_directories(AutoConvertProperties
  PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}")

target_link_libraries(AutoConvertProperties
  PRIVATE
    ParaView::pqApplicationComponents
    ParaView::RemotingServerManager
)

install(TARGETS AutoConvertProperties
  RUNTIME DESTINATION lib/paraview
  LIBRARY DESTINATION lib/paraview
  ARCHIVE DESTINATION lib/paraview
)

