#include <benchmark/benchmark.h>

class DomainBenchmark
{
public:
    int nx = 50;
    int ny = 50;
    int nz = 50;
    int nyTimesnz;

    DomainBenchmark() : nyTimesnz(ny * nz) {}

    int getGlobalIndexOriginal(const int &x, const int &y, const int &z)
    {
        return z + nz * (y + ny * x);
    }

    int getGlobalIndexNoRef(int x, int y, int z)
    {
        return z + nz * (y + ny * x);
    }

    int getGlobalIndexPrecalc(int x, int y, int z)
    {
        return z + nz * y + nyTimesnz * x;
    }
};

static void BM_GlobalIndexOriginal(benchmark::State &state)
{
    DomainBenchmark domain;
    for (auto _ : state)
    {
        for (int i = 0; i < domain.nx; i++)
        {
            for (int j = 0; j < domain.ny; j++)
            {
                for (int k = 0; k < domain.nz; k++)
                {
                    benchmark::DoNotOptimize(domain.getGlobalIndexOriginal(i, j, k));
                }
            }
        }
    }
}
BENCHMARK(BM_GlobalIndexOriginal);

static void BM_GlobalIndexNoRef(benchmark::State &state)
{
    DomainBenchmark domain;
    for (auto _ : state)
    {
        for (int i = 0; i < domain.nx; i++)
        {
            for (int j = 0; j < domain.ny; j++)
            {
                for (int k = 0; k < domain.nz; k++)
                {
                    benchmark::DoNotOptimize(domain.getGlobalIndexNoRef(i, j, k));
                }
            }
        }
    }
}
BENCHMARK(BM_GlobalIndexNoRef);

static void BM_GlobalIndexPrecalc(benchmark::State &state)
{
    DomainBenchmark domain;
    for (auto _ : state)
    {
        for (int i = 0; i < domain.nx; i++)
        {
            for (int j = 0; j < domain.ny; j++)
            {
                for (int k = 0; k < domain.nz; k++)
                {
                    benchmark::DoNotOptimize(domain.getGlobalIndexPrecalc(i, j, k));
                }
            }
        }
    }
}
BENCHMARK(BM_GlobalIndexPrecalc);

BENCHMARK_MAIN();