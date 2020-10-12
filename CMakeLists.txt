cmake_minimum_required(VERSION 3.16)
project(GeodesicPath)

set(CMAKE_CXX_STANDARD 17)

set(INCLUDE /usr/local/include)
set(SURFACE_MESH /usr/local/include/SurfaceMesh)
set(SURFACE_MESH_LIB /usr/local/lib/libSurfaceMesh.dylib)

include_directories(${SURFACE_MESH} ${INCLUDE})

# Eigen
find_package(Eigen3 3.3.0)
if (Eigen3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
else ()
    message(can't find eigen3)
    return()
endif ()

aux_source_directory(. DIR_SRC)

add_executable(${PROJECT_NAME} ${DIR_SRC})

target_link_libraries(${PROJECT_NAME} ${SURFACE_MESH_LIB})