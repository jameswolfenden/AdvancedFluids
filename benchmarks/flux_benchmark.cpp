#include <benchmark/benchmark.h>
#include "types/Flux.hpp"

// Original version
struct FluxOriginal
{
    double f1, f2, f3, f4, f5;
    FluxOriginal() : f1(0), f2(0), f3(0), f4(0), f5(0) {}
    void updateFromPrimatives(double rho, double u, double v, double w, double p)
    {
        f1 = rho;
        f2 = rho * u;
        f3 = rho * v;
        f4 = rho * w;
        f5 = 0.5 * rho * (u * u + v * v + w * w) + p / (1.4 - 1);
    }
};

// 32-byte aligned version
struct alignas(32) Flux32
{
    double f1, f2, f3, f4, f5;
    Flux32() : f1(0), f2(0), f3(0), f4(0), f5(0) {}
    void updateFromPrimatives(double rho, double u, double v, double w, double p)
    {
        f1 = rho;
        f2 = rho * u;
        f3 = rho * v;
        f4 = rho * w;
        f5 = 0.5 * rho * (u * u + v * v + w * w) + p / (1.4 - 1);
    }
};

// 64-byte aligned version
struct alignas(64) Flux64
{
    double f1, f2, f3, f4, f5;
    Flux64() : f1(0), f2(0), f3(0), f4(0), f5(0) {}
    void updateFromPrimatives(double rho, double u, double v, double w, double p)
    {
        f1 = rho;
        f2 = rho * u;
        f3 = rho * v;
        f4 = rho * w;
        f5 = 0.5 * rho * (u * u + v * v + w * w) + p / (1.4 - 1);
    }
};

static void BM_FluxOriginal(benchmark::State &state)
{
    std::vector<FluxOriginal> fluxes(1000);
    for (auto _ : state)
    {
        for (auto &f : fluxes)
        {
            f.updateFromPrimatives(1, 2, 3, 4, 5);
            benchmark::DoNotOptimize(f);
        }
    }
}
BENCHMARK(BM_FluxOriginal);

static void BM_Flux32(benchmark::State &state)
{
    std::vector<Flux32> fluxes(1000);
    for (auto _ : state)
    {
        for (auto &f : fluxes)
        {
            f.updateFromPrimatives(1, 2, 3, 4, 5);
            benchmark::DoNotOptimize(f);
        }
    }
}
BENCHMARK(BM_Flux32);

static void BM_Flux64(benchmark::State &state)
{
    std::vector<Flux64> fluxes(1000);
    for (auto _ : state)
    {
        for (auto &f : fluxes)
        {
            f.updateFromPrimatives(1, 2, 3, 4, 5);
            benchmark::DoNotOptimize(f);
        }
    }
}
BENCHMARK(BM_Flux64);

BENCHMARK_MAIN();