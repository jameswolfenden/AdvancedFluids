#find_package(benchmark REQUIRED)

add_executable(index_benchmark getGlobalIndex_benchmark.cpp)
target_link_libraries(index_benchmark 
    PRIVATE 
    benchmark::benchmark 
    fluids_lib    # Your main library target
)
