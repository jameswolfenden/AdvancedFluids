#include <benchmark/benchmark.h>
#include <omp.h>
#include "solvers/DomainEulerSolver.hpp"
#include "types/State.hpp"

namespace fluid
{

    class DomainEulerSolverBenchmark
    {
    public:
        std::vector<Domain> domains;

        DomainEulerSolverBenchmark() : solver(0.7)
        {
            omp_set_num_threads(6);
            setupBenchmarkDomains(domains);
        }

        // Methods to directly access private functions for benchmarking
        bool benchmarkTimeStep()
        {
            solver.fetchTimeStep(domains);
            return true;
        }

        bool benchmarkGhostCells()
        {
            solver.updateGhostCells(domains);
            return true;
        }

        bool benchmarkXFaces(Domain &domain)
        {
            return solver.xFaces(domain);
        }

        bool benchmarkYFaces(Domain &domain)
        {
            return solver.yFaces(domain);
        }

        bool benchmarkZFaces(Domain &domain)
        {
            return solver.zFaces(domain);
        }

    private:
        static void setupBenchmarkDomains(std::vector<Domain> &domains)
        {
            // Create simple test case with 4 domains
            domains.resize(4);

            // Initial conditions, atmospheric conditions
            State initial(1.0, 0.0, 0.0, 0.0, 101325.0);

            // Setup domains with 50x50x50 cells each
            double cellDensity = 50;
            double boxSize = 1.0;

            for (int i = 0; i < 4; i++)
            {
                domains[i].setup(i, boxSize, boxSize, boxSize, cellDensity, initial);
            }

            // Connect domains in 2x2 configuration
            // Domain 0 connects to 1 on right and 2 below
            domains[0].sides[0] = &domains[1]; // +x side
            domains[0].sides[2] = &domains[2]; // +y side

            // Domain 1 connects to 0 on left and 3 below
            domains[1].sides[1] = &domains[0]; // -x side
            domains[1].sides[2] = &domains[3]; // +y side

            // Domain 2 connects to 0 above and 3 on right
            domains[2].sides[3] = &domains[0]; // -y side
            domains[2].sides[0] = &domains[3]; // +x side

            // Domain 3 connects to 1 above and 2 on left
            domains[3].sides[3] = &domains[1]; // -y side
            domains[3].sides[1] = &domains[2]; // -x side
        }
        DomainEulerSolver solver;
    };

} // namespace fluid

// Benchmark definitions
static void BM_TimeStepCalculation(benchmark::State &state)
{
    fluid::DomainEulerSolverBenchmark bench;

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bench.benchmarkTimeStep());
    }
}
BENCHMARK(BM_TimeStepCalculation);

static void BM_GhostCellUpdate(benchmark::State &state)
{
    fluid::DomainEulerSolverBenchmark bench;

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bench.benchmarkGhostCells());
    }
}
BENCHMARK(BM_GhostCellUpdate);

static void BM_XFacesCalculation(benchmark::State &state)
{
    fluid::DomainEulerSolverBenchmark bench;

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bench.benchmarkXFaces(bench.domains[0]));
    }
}
BENCHMARK(BM_XFacesCalculation);

static void BM_YFacesCalculation(benchmark::State &state)
{
    fluid::DomainEulerSolverBenchmark bench;

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bench.benchmarkYFaces(bench.domains[0]));
    }
}
BENCHMARK(BM_YFacesCalculation);

static void BM_ZFacesCalculation(benchmark::State &state)
{
    fluid::DomainEulerSolverBenchmark bench;

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bench.benchmarkZFaces(bench.domains[0]));
    }
}
BENCHMARK(BM_ZFacesCalculation);

BENCHMARK_MAIN();