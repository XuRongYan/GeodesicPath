set(GOOGLE_TEST_ROOT external/googletest)

include_directories(
        ${PROJECT_SOURCE_DIR}/${GOOGLE_TEST_ROOT}
        ${PROJECT_SOURCE_DIR}/${GOOGLE_TEST_ROOT}/include
)

message(${GOOGLE_TEST_ROOT})

set(GOOGLE_TEST_SOURCES
        ${PROJECT_SOURCE_DIR}/${GOOGLE_TEST_ROOT}/src/gtest-all.cc
        ${PROJECT_SOURCE_DIR}/${GOOGLE_TEST_ROOT}/src/gtest_main.cc
        )

add_library(googletest ${GOOGLE_TEST_SOURCES})

# set output path so that coverage testing works
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

# source files
file(GLOB SOURCES ./*.cpp)

aux_source_directory(../ SRC)

# build test runner
add_executable(gtest_runner ${SOURCES} ${SRC})

target_link_libraries(gtest_runner
        ${OPENMP}
        ${SURFACE_MESH_LIB}
        )

if(UNIX OR APPLE)
    target_link_libraries(gtest_runner googletest pthread)
elseif(WIN32)
    target_link_libraries(gtest_runner googletest)
endif()

add_dependencies(gtest_runner googletest)

# add runner as test
add_test(gtest_runner ${CMAKE_CURRENT_BINARY_DIR}/gtest_runner)