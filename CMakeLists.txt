cmake_minimum_required(VERSION 3.16)
project(my_tet_arap)

set(CMAKE_CXX_STANDARD 11)

set(VTK_COMMON_CORE /usr/local/lib/libvtkCommonCore-8.2.dylib)
set(VTK_COMMON_DATA_MODEL /usr/local/lib/libvtkCommonDataModel-8.2.dylib)
set(VTK_IO_CORE /usr/local/lib/libvtkIOCore-8.2.dylib)
set(VTK_IO_LEGACY /usr/local/lib/libvtkIOLegacy-8.2.dylib)
set(VTK_IO_XML /usr/local/lib/libvtkIOXML-8.2.dylib)

include_directories(/usr/local/include)
include_directories(/usr/local/include/eigen3)
include_directories(/usr/local/include/SurfaceMesh)
include_directories(/usr/local/include/vtk-8.2)

add_compile_definitions(CINOLIB_USES_VTK)

aux_source_directory(. DIRS_SRCS)

add_executable(my_tet_arap ${DIRS_SRCS})

target_link_libraries(my_tet_arap ${SURFACE_MESH} ${VTK_COMMON_CORE} ${VTK_COMMON_DATA_MODEL} ${VTK_IO_CORE} ${VTK_IO_LEGACY} ${VTK_IO_XML})