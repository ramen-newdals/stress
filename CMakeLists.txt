cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

project(stress C CXX Fortran)

#===========HDF5 CONFIGURATION =======================
# set HDF5 source dir
set(HDF5_DIR "/home/ramen_newdals/CMake-hdf5-1.12.2/build/HDF5-1.12.2-Linux/HDF_Group/HDF5/1.12.2/cmake")

set (FIND_HDF_COMPONENTS C CXX shared)

find_package (HDF5 NAMES "hdf5" COMPONENTS ${FIND_HDF_COMPONENTS})
if (HDF5_FOUND)
  if (HDF5_shared_C_FOUND)
    set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_SHARED_LIBRARY})
  endif ()
  if (HDF5_shared_CXX_FOUND)
    set (LINK_LIBS ${LINK_LIBS} ${HDF5_CXX_SHARED_LIBRARY})
  endif ()
else ()
  message (FATAL_ERROR " HDF5 is Required")
endif ()
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})


#=========VTK CONFIGURATION ========================
find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonDataModel
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

set (LINK_LIBS ${LINK_LIBS} ${VTK_LIBRARIES})

#=========nl CONFIGURATION==========================
set (OPENNL_SRCS
    nl_single_file.c
    )

include_directories(${OPENNL_SOURCE_DIR})

add_library (nl STATIC ${OPENNL_SRCS})

if(NOT WIN32)
  set_target_properties(nl PROPERTIES COMPILE_FLAGS -fPIC)
endif()

install(TARGETS nl
  EXPORT VMTK-Targets
  LIBRARY DESTINATION ${VTK_VMTK_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${VTK_VMTK_INSTALL_LIB_DIR} COMPONENT Development
  RUNTIME DESTINATION ${VTK_VMTK_INSTALL_BIN_DIR} COMPONENT RuntimeExecutables
)

install(FILES nl.h
  DESTINATION ${VTK_VMTK_INSTALL_INCLUDE_DIR}
  COMPONENT Development)


#==========Linking and generation of Executable ===========
message("we are in the beam")
message(${HDF5_INCLUDE_DIR})

add_executable(readMesh readMesh.cpp gq.cpp)

target_link_libraries (readMesh PRIVATE ${LINK_LIBS})