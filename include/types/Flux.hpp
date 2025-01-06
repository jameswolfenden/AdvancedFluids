#ifndef FLUX_HPP
#define FLUX_HPP

#include "constants.hpp"

namespace fluid
{

    struct Flux
    {
        double f1, f2, f3, f4, f5;

        Flux();
        Flux(double f1, double f2, double f3, double f4, double f5);
        void updateFromPrimatives(double rho, double u, double v, double w, double p);
    };

} // namespace fluid

#endif // FLUX_HPP