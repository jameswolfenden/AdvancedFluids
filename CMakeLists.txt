cmake_minimum_required(VERSION 3.5.0)
project(fluids VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

include_directories(${INCLUDE_DIR})

# Collect all source files except main.cpp
file(GLOB_RECURSE SOURCES 
    "${SRC_DIR}/*.cpp"
)
list(FILTER SOURCES EXCLUDE REGEX ".*/main\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/main_old\\.cpp$")

# Print found sources
message(STATUS "Found source files:")
foreach(SOURCE ${SOURCES})
    message(STATUS "  ${SOURCE}")
endforeach()

# Create library from sources
add_library(fluids_lib STATIC ${SOURCES})

# Create main executable
add_executable(fluids ${SRC_DIR}/main.cpp)

# Find and link HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)
if (HDF5_FOUND)
    # Print HDF5 information for debugging
    message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
    message(STATUS "HDF5_CXX_INCLUDE_DIRS: ${HDF5_CXX_INCLUDE_DIRS}")
    message(STATUS "HDF5_CXX_LIBRARIES: ${HDF5_CXX_LIBRARIES}")

    # Link to library
    target_include_directories(fluids_lib PUBLIC ${HDF5_INCLUDE_DIRS} ${HDF5_CXX_INCLUDE_DIRS})
    target_link_libraries(fluids_lib PUBLIC ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})
endif()

# Find and link OpenMP
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    target_link_libraries(fluids_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Link executable to our library (dependencies will propagate)
target_link_libraries(fluids PRIVATE fluids_lib)

# Add compiler warnings
if(MSVC)
    target_compile_options(fluids_lib PRIVATE /W4)
    target_compile_options(fluids PRIVATE /W4)
else()
    target_compile_options(fluids_lib PRIVATE -Wall -Wextra)
    target_compile_options(fluids PRIVATE -Wall -Wextra)
endif()

# Optional: Add tests if they exist
if(EXISTS "${PROJECT_SOURCE_DIR}/tests/CMakeLists.txt")
    enable_testing()
    add_subdirectory(tests)
endif()

# Print status information
message(STATUS "HDF5 libs: ${HDF5_LIBRARIES}")
message(STATUS "HDF5 includes: ${HDF5_INCLUDE_DIRS}")