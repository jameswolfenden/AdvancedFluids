find_package(benchmark REQUIRED)

add_executable(index_benchmark getGlobalIndex_benchmark.cpp)
target_link_libraries(index_benchmark 
    PRIVATE 
    benchmark::benchmark 
    fluids_lib
)

add_executable(euler_benchmark domainEulerSolver_benchmark.cpp)
target_link_libraries(euler_benchmark 
    PRIVATE 
    benchmark::benchmark 
    fluids_lib
)

add_executable(flux_benchmark flux_benchmark.cpp)
target_link_libraries(flux_benchmark 
    PRIVATE 
    benchmark::benchmark 
    fluids_lib
)

add_executable(euler_benchmark_parallel parallelDomainEulerSolver_benchmark.cpp)
target_link_libraries(euler_benchmark_parallel 
    PRIVATE 
    benchmark::benchmark 
    fluids_lib
)
