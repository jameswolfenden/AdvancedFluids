# Create test executable
add_executable(fluids_tests
    solvers/RiemannSolverTest.cpp
)

# Link against Google Test and your library
target_link_libraries(fluids_tests
    PRIVATE
    GTest::gtest_main
    fluids_lib
)

# Enable test discovery
include(GoogleTest)
gtest_discover_tests(fluids_tests)

#add_test(NAME fluids_tests COMMAND fluids_tests)
