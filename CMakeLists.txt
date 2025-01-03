cmake_minimum_required(VERSION 3.5.0)
project(fluids VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

include_directories(${INCLUDE_DIR})

add_executable(fluids ${SRC_DIR}/main.cpp)

find_package(HDF5 REQUIRED)
if (HDF5_FOUND)
    target_include_directories(fluids PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(fluids PRIVATE ${HDF5_LIBRARIES})
endif()

find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    target_link_libraries(fluids PRIVATE OpenMP::OpenMP_CXX)
endif()
