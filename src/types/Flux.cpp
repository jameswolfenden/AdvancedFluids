#include "types/Flux.hpp"

namespace fluid
{

    Flux::Flux()
        : f1(0.0), f2(0.0), f3(0.0), f4(0.0), f5(0.0) {}

    Flux::Flux(double f1, double f2, double f3, double f4, double f5)
        : f1(f1), f2(f2), f3(f3), f4(f4), f5(f5) {}

    void Flux::updateFromPrimatives(double rho, double u, double v, double w, double p)
    {
        f1 = rho * u;
        f2 = rho * u * u + p;
        f3 = rho * v * u;
        f4 = rho * w * u;
        f5 = u * (rho * (0.5 * (u * u + v * v + w * w) + p / ((G - 1) * rho) + p));
    }

    void Flux::updateFromPrimatives(const StateRef &state)
    {
        updateFromPrimatives(state.rho, state.u, state.v, state.w, state.p);
    }

} // namespace fluid