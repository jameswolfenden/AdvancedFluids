#include <benchmark/benchmark.h>
#include <omp.h>
#include "types/ConservedState.hpp"
#include "types/State.hpp"
#include "domain/Domain.hpp"
#include "solvers/RiemannSolver.hpp"
#include "domain/DomainPositioner.hpp"

namespace fluid
{
    class DomainEulerSolver_original
    {
    public:
        RiemannSolver rs;
        double minT = 0.001; // Fixed timestep for testing
        DomainEulerSolver_original()
        {
            rs = RiemannSolver();
        }

        bool updateX(std::vector<Domain> &domains)
        {
            std::atomic<bool> errorFlag(false);
#pragma omp parallel for
            for (auto &domain : domains)
            {
                if (!xFaces(domain))
                    errorFlag.store(true);
                if (!xBoxes(domain))
                    errorFlag.store(true);
            }
            return !errorFlag.load();
        }
        bool xFaces(Domain &domain)
        {
            std::atomic<bool> errorFlag(false);
#pragma omp parallel for collapse(3)
            for (int i = 0; i < domain.nx; i++)
            {
                for (int j = 0; j < domain.ny; j++)
                {
                    for (int k = 0; k < domain.nzFaces; k++)
                    {
                        if (errorFlag.load())
                            continue;

                        int index = k + domain.nz * (j + domain.ny * i);

                        StateRef left(domain.rho_[index], domain.w_[index], domain.v_[index],
                                       domain.u_[index], domain.p_[index], domain.a_[index]);

                        index++;

                        StateRef right(domain.rho_[index], domain.w_[index], domain.v_[index],
                                        domain.u_[index], domain.p_[index], domain.a_[index]);

                        if (!rs.findStar(left, right, domain.xfAt(i, j, k)))
                        {

                            errorFlag.store(true);
                        }
                    }
                }
            }
            return !errorFlag.load();
        }
        bool xBoxes(Domain &domain)
        {
            std::atomic<bool> errorFlag(false);
#pragma omp parallel for collapse(3)
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int j = 1; j < domain.ny - 1; j++)
                {
                    for (int k = 1; k < domain.nz - 1; k++)
                    {
                        if (errorFlag.load())
                            continue;

                        ConservedState state(domain.rho(i, j, k), domain.u(i, j, k),
                                             domain.v(i, j, k), domain.w(i, j, k),
                                             domain.p(i, j, k), domain.a(i, j, k));

                        double u1 = state.u1() + minT * (domain.xfAt(i - 1, j, k).f1 - domain.xfAt(i, j, k).f1) / domain.boxDims;
                        double u2 = state.u2() + minT * (domain.xfAt(i - 1, j, k).f2 - domain.xfAt(i, j, k).f2) / domain.boxDims;
                        double u3 = state.u3() + minT * (domain.xfAt(i - 1, j, k).f3 - domain.xfAt(i, j, k).f3) / domain.boxDims;
                        double u4 = state.u4() + minT * (domain.xfAt(i - 1, j, k).f4 - domain.xfAt(i, j, k).f4) / domain.boxDims;
                        double u5 = state.u5() + minT * (domain.xfAt(i - 1, j, k).f5 - domain.xfAt(i, j, k).f5) / domain.boxDims;

                        if (!state.updateFromConservatives(u1, u2, u3, u4, u5))
                        {
                            errorFlag.store(true);
                        }
                    }
                }
            }
            return !errorFlag.load();
        }
    };

    class DomainEulerSolver_domainLevel
    {
    public:
        RiemannSolver rs;
        double minT = 0.001; // Fixed timestep for testing
        DomainEulerSolver_domainLevel()
        {
            rs = RiemannSolver();
        }

        bool updateX(std::vector<Domain> &domains)
        {
            std::atomic<bool> errorFlag(false);
#pragma omp parallel for
            for (auto &domain : domains)
            {
                if (!xFaces(domain))
                    errorFlag.store(true);
                if (!xBoxes(domain))
                    errorFlag.store(true);
            }
            return !errorFlag.load();
        }
        bool xFaces(Domain &domain)
        {
            std::atomic<bool> errorFlag(false);
            for (int i = 0; i < domain.nx; i++)
            {
                for (int j = 0; j < domain.ny; j++)
                {
                    for (int k = 0; k < domain.nzFaces; k++)
                    {
                        if (errorFlag.load())
                            continue;

                        int index = k + domain.nz * (j + domain.ny * i);

                        StateRef left(domain.rho_[index], domain.w_[index], domain.v_[index],
                                       domain.u_[index], domain.p_[index], domain.a_[index]);

                        index++;

                        StateRef right(domain.rho_[index], domain.w_[index], domain.v_[index],
                                        domain.u_[index], domain.p_[index], domain.a_[index]);

                        if (!rs.findStar(left, right, domain.xfAt(i, j, k)))
                        {

                            errorFlag.store(true);
                        }
                    }
                }
            }
            return !errorFlag.load();
        }
        bool xBoxes(Domain &domain)
        {
            std::atomic<bool> errorFlag(false);
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int j = 1; j < domain.ny - 1; j++)
                {
                    for (int k = 1; k < domain.nz - 1; k++)
                    {
                        if (errorFlag.load())
                            continue;

                        ConservedState state(domain.rho(i, j, k), domain.u(i, j, k),
                                             domain.v(i, j, k), domain.w(i, j, k),
                                             domain.p(i, j, k), domain.a(i, j, k));

                        double u1 = state.u1() + minT * (domain.xfAt(i - 1, j, k).f1 - domain.xfAt(i, j, k).f1) / domain.boxDims;
                        double u2 = state.u2() + minT * (domain.xfAt(i - 1, j, k).f2 - domain.xfAt(i, j, k).f2) / domain.boxDims;
                        double u3 = state.u3() + minT * (domain.xfAt(i - 1, j, k).f3 - domain.xfAt(i, j, k).f3) / domain.boxDims;
                        double u4 = state.u4() + minT * (domain.xfAt(i - 1, j, k).f4 - domain.xfAt(i, j, k).f4) / domain.boxDims;
                        double u5 = state.u5() + minT * (domain.xfAt(i - 1, j, k).f5 - domain.xfAt(i, j, k).f5) / domain.boxDims;

                        if (!state.updateFromConservatives(u1, u2, u3, u4, u5))
                        {
                            errorFlag.store(true);
                        }
                    }
                }
            }
            return !errorFlag.load();
        }
    };

    class DomainEulerSolver_loopLevel
    {
    public:
        RiemannSolver rs;
        double minT = 0.001; // Fixed timestep for testing
        DomainEulerSolver_loopLevel()
        {
            rs = RiemannSolver();
        }

        bool updateX(std::vector<Domain> &domains)
        {
            std::atomic<bool> errorFlag(false);
            for (auto &domain : domains)
            {
                if (!xFaces(domain))
                    errorFlag.store(true);
                if (!xBoxes(domain))
                    errorFlag.store(true);
            }
            return !errorFlag.load();
        }
        bool xFaces(Domain &domain)
        {
            std::atomic<bool> errorFlag(false);
#pragma omp parallel for collapse(3)
            for (int i = 0; i < domain.nx; i++)
            {
                for (int j = 0; j < domain.ny; j++)
                {
                    for (int k = 0; k < domain.nzFaces; k++)
                    {
                        if (errorFlag.load())
                            continue;

                        int index = k + domain.nz * (j + domain.ny * i);

                        StateRef left(domain.rho_[index], domain.w_[index], domain.v_[index],
                                       domain.u_[index], domain.p_[index], domain.a_[index]);

                        index++;

                        StateRef right(domain.rho_[index], domain.w_[index], domain.v_[index],
                                        domain.u_[index], domain.p_[index], domain.a_[index]);

                        if (!rs.findStar(left, right, domain.xfAt(i, j, k)))
                        {

                            errorFlag.store(true);
                        }
                    }
                }
            }
            return !errorFlag.load();
        }
        bool xBoxes(Domain &domain)
        {
            std::atomic<bool> errorFlag(false);
#pragma omp parallel for collapse(3)
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int j = 1; j < domain.ny - 1; j++)
                {
                    for (int k = 1; k < domain.nz - 1; k++)
                    {
                        if (errorFlag.load())
                            continue;

                        ConservedState state(domain.rho(i, j, k), domain.u(i, j, k),
                                             domain.v(i, j, k), domain.w(i, j, k),
                                             domain.p(i, j, k), domain.a(i, j, k));

                        double u1 = state.u1() + minT * (domain.xfAt(i - 1, j, k).f1 - domain.xfAt(i, j, k).f1) / domain.boxDims;
                        double u2 = state.u2() + minT * (domain.xfAt(i - 1, j, k).f2 - domain.xfAt(i, j, k).f2) / domain.boxDims;
                        double u3 = state.u3() + minT * (domain.xfAt(i - 1, j, k).f3 - domain.xfAt(i, j, k).f3) / domain.boxDims;
                        double u4 = state.u4() + minT * (domain.xfAt(i - 1, j, k).f4 - domain.xfAt(i, j, k).f4) / domain.boxDims;
                        double u5 = state.u5() + minT * (domain.xfAt(i - 1, j, k).f5 - domain.xfAt(i, j, k).f5) / domain.boxDims;

                        if (!state.updateFromConservatives(u1, u2, u3, u4, u5))
                        {
                            errorFlag.store(true);
                        }
                    }
                }
            }
            return !errorFlag.load();
        }
    };

    class DomainEulerSolver_none
    {
    public:
        RiemannSolver rs;
        double minT = 0.001; // Fixed timestep for testing
        DomainEulerSolver_none()
        {
            rs = RiemannSolver();
        }

        bool updateX(std::vector<Domain> &domains)
        {
            std::atomic<bool> errorFlag(false);
            for (auto &domain : domains)
            {
                if (!xFaces(domain))
                    errorFlag.store(true);
                if (!xBoxes(domain))
                    errorFlag.store(true);
            }
            return !errorFlag.load();
        }
        bool xFaces(Domain &domain)
        {
            std::atomic<bool> errorFlag(false);
            for (int i = 0; i < domain.nx; i++)
            {
                for (int j = 0; j < domain.ny; j++)
                {
                    for (int k = 0; k < domain.nzFaces; k++)
                    {
                        if (errorFlag.load())
                            continue;

                        int index = k + domain.nz * (j + domain.ny * i);

                        StateRef left(domain.rho_[index], domain.w_[index], domain.v_[index],
                                       domain.u_[index], domain.p_[index], domain.a_[index]);

                        index++;

                        StateRef right(domain.rho_[index], domain.w_[index], domain.v_[index],
                                        domain.u_[index], domain.p_[index], domain.a_[index]);

                        if (!rs.findStar(left, right, domain.xfAt(i, j, k)))
                        {

                            errorFlag.store(true);
                        }
                    }
                }
            }
            return !errorFlag.load();
        }
        bool xBoxes(Domain &domain)
        {
            std::atomic<bool> errorFlag(false);
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int j = 1; j < domain.ny - 1; j++)
                {
                    for (int k = 1; k < domain.nz - 1; k++)
                    {
                        if (errorFlag.load())
                            continue;

                        ConservedState state(domain.rho(i, j, k), domain.u(i, j, k),
                                             domain.v(i, j, k), domain.w(i, j, k),
                                             domain.p(i, j, k), domain.a(i, j, k));

                        double u1 = state.u1() + minT * (domain.xfAt(i - 1, j, k).f1 - domain.xfAt(i, j, k).f1) / domain.boxDims;
                        double u2 = state.u2() + minT * (domain.xfAt(i - 1, j, k).f2 - domain.xfAt(i, j, k).f2) / domain.boxDims;
                        double u3 = state.u3() + minT * (domain.xfAt(i - 1, j, k).f3 - domain.xfAt(i, j, k).f3) / domain.boxDims;
                        double u4 = state.u4() + minT * (domain.xfAt(i - 1, j, k).f4 - domain.xfAt(i, j, k).f4) / domain.boxDims;
                        double u5 = state.u5() + minT * (domain.xfAt(i - 1, j, k).f5 - domain.xfAt(i, j, k).f5) / domain.boxDims;

                        if (!state.updateFromConservatives(u1, u2, u3, u4, u5))
                        {
                            errorFlag.store(true);
                        }
                    }
                }
            }
            return !errorFlag.load();
        }
    };

    class ParallelizationBenchmark
    {
    public:
        std::vector<Domain> domains;
        DomainEulerSolver_original original;
        DomainEulerSolver_domainLevel domainLevel;
        DomainEulerSolver_loopLevel loopLevel;
        DomainEulerSolver_none none;

        ParallelizationBenchmark(int numDomains, int cellsPerDim)
        {
            setupBenchmarkDomains(numDomains, cellsPerDim);
        }

        bool runOriginal()
        {
            original.updateX(domains);
            return true;
        }

        bool runDomainLevel()
        {
            domainLevel.updateX(domains);
            return true;
        }

        bool runLoopLevel()
        {
            loopLevel.updateX(domains);
            return true;
        }

        bool runNone()
        {
            none.updateX(domains);
            return true;
        }

    private:
        void setupBenchmarkDomains(int numDomains, int cellsPerDim)
        {
            domains.resize(numDomains);

            // Initial conditions - standard atmospheric conditions
            State initial(1.0, 0.0, 0.0, 0.0, 101325.0);

            // Set up domains with specified size
            for (int i = 0; i < numDomains; i++)
            {
                domains[i].setup(i, 1.0, 1.0, 1.0, cellsPerDim, initial);
            }

            // Connect domains in the x-direction
            for (int i = 0; i < numDomains - 1; i++)
            {
                domains[i].sides[0] = &domains[i + 1];
                domains[i + 1].sides[1] = &domains[i];
            }

        }
    };
} // namespace fluid

// Benchmark with different domain counts and sizes
static void BM_Original(benchmark::State &state)
{
    const int numDomains = state.range(0);  // Number of domains
    const int cellsPerDim = state.range(1); // Cells per dimension
    const int numThreads = state.range(2);  // Number of threads

    fluid::ParallelizationBenchmark bench(numDomains, cellsPerDim);
    omp_set_num_threads(numThreads);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bench.runOriginal());
    }

    // Custom counter for problem size
    state.counters["Domains"] = numDomains;
    state.counters["Cells_Per_Dim"] = cellsPerDim;
    state.counters["Total_Cells"] = numDomains * cellsPerDim * cellsPerDim * cellsPerDim;
    state.counters["Threads"] = numThreads;
}

static void BM_DomainLevel(benchmark::State &state)
{
    const int numDomains = state.range(0);
    const int cellsPerDim = state.range(1);
    const int numThreads = state.range(2);

    fluid::ParallelizationBenchmark bench(numDomains, cellsPerDim);
    omp_set_num_threads(numThreads);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bench.runDomainLevel());
    }

    state.counters["Domains"] = numDomains;
    state.counters["Cells_Per_Dim"] = cellsPerDim;
    state.counters["Total_Cells"] = numDomains * cellsPerDim * cellsPerDim * cellsPerDim;
    state.counters["Threads"] = numThreads;
}

static void BM_LoopLevel(benchmark::State &state)
{
    const int numDomains = state.range(0);
    const int cellsPerDim = state.range(1);
    const int numThreads = state.range(2);

    fluid::ParallelizationBenchmark bench(numDomains, cellsPerDim);
    omp_set_num_threads(numThreads);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bench.runLoopLevel());
    }

    state.counters["Domains"] = numDomains;
    state.counters["Cells_Per_Dim"] = cellsPerDim;
    state.counters["Total_Cells"] = numDomains * cellsPerDim * cellsPerDim * cellsPerDim;
    state.counters["Threads"] = numThreads;
}

static void BM_None(benchmark::State &state)
{
    const int numDomains = state.range(0);
    const int cellsPerDim = state.range(1);
    const int numThreads = state.range(2);

    fluid::ParallelizationBenchmark bench(numDomains, cellsPerDim);
    omp_set_num_threads(numThreads);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bench.runNone());
    }

    state.counters["Domains"] = numDomains;
    state.counters["Cells_Per_Dim"] = cellsPerDim;
    state.counters["Total_Cells"] = numDomains * cellsPerDim * cellsPerDim * cellsPerDim;
    state.counters["Threads"] = numThreads;
}

// Register benchmarks with various configurations
// Args: (num_domains, cells_per_dim, num_threads)

// Test with few large domains
BENCHMARK(BM_Original)->Args({2, 100, 1})->Args({2, 100, 2})->Args({2, 100, 4})->Args({2, 100, 8})->Args({4, 100, 1})->Args({4, 100, 2})->Args({4, 100, 4})->Args({4, 100, 8})->Unit(benchmark::kMillisecond);

BENCHMARK(BM_DomainLevel)->Args({2, 100, 1})->Args({2, 100, 2})->Args({2, 100, 4})->Args({2, 100, 8})->Args({4, 100, 1})->Args({4, 100, 2})->Args({4, 100, 4})->Args({4, 100, 8})->Unit(benchmark::kMillisecond);

BENCHMARK(BM_LoopLevel)->Args({2, 100, 1})->Args({2, 100, 2})->Args({2, 100, 4})->Args({2, 100, 8})->Args({4, 100, 1})->Args({4, 100, 2})->Args({4, 100, 4})->Args({4, 100, 8})->Unit(benchmark::kMillisecond);

BENCHMARK(BM_None)->Args({2, 100, 1})->Args({2, 100, 2})->Args({2, 100, 4})->Args({2, 100, 8})->Args({4, 100, 1})->Args({4, 100, 2})->Args({4, 100, 4})->Args({4, 100, 8})->Unit(benchmark::kMillisecond);

// Test with many small domains
BENCHMARK(BM_Original)->Args({8, 50, 1})->Args({8, 50, 2})->Args({8, 50, 4})->Args({8, 50, 8})->Args({16, 50, 1})->Args({16, 50, 2})->Args({16, 50, 4})->Args({16, 50, 8})->Unit(benchmark::kMillisecond);

BENCHMARK(BM_DomainLevel)->Args({8, 50, 1})->Args({8, 50, 2})->Args({8, 50, 4})->Args({8, 50, 8})->Args({16, 50, 1})->Args({16, 50, 2})->Args({16, 50, 4})->Args({16, 50, 8})->Unit(benchmark::kMillisecond);

BENCHMARK(BM_LoopLevel)->Args({8, 50, 1})->Args({8, 50, 2})->Args({8, 50, 4})->Args({8, 50, 8})->Args({16, 50, 1})->Args({16, 50, 2})->Args({16, 50, 4})->Args({16, 50, 8})->Unit(benchmark::kMillisecond);

BENCHMARK(BM_None)->Args({8, 50, 1})->Args({8, 50, 2})->Args({8, 50, 4})->Args({8, 50, 8})->Args({16, 50, 1})->Args({16, 50, 2})->Args({16, 50, 4})->Args({16, 50, 8})->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
