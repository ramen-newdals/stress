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

message("we are in the beam")
message(${HDF5_INCLUDE_DIR})

add_executable(readMesh readMesh.cpp gq.cpp)

target_link_libraries (readMesh PRIVATE ${LINK_LIBS})