#include "types/ConservedState.hpp"
#include "constants.hpp"
#include <cmath>
#include <iostream>

namespace fluid
{

    ConservedState::ConservedState(double &rho, double &u, double &v, double &w, double &p, double &a)
        : rho(rho), u(u), v(v), w(w), p(p), a(a) {}

    double ConservedState::aCalc()
    {
        if (rho == 0.0)
        {
            a = 0.0;
        }
        else
        {
            a = sqrt((G * p) / rho);
        }
        return a;
    }

    void ConservedState::fixVacuum()
    {
        if (rho <= 0.0)
        {
            p = 0.0;
        }
        if (p <= 0.0)
        {
            rho = 0.0;
            p = 0.0;
        }
        aCalc();
    }

    double ConservedState::u1() const
    {
        return rho;
    }

    double ConservedState::u2() const
    {
        return rho * u;
    }

    double ConservedState::u3() const
    {
        return rho * v;
    }

    double ConservedState::u4() const
    {
        return rho * w;
    }

    double ConservedState::u5() const
    {
        if (rho == 0.0)
            return 0.0;
        return rho * 0.5 * (u * u + v * v + w * w) + p / (G - 1);
    }

    bool ConservedState::updateFromConservatives(double u1, double u2, double u3, double u4, double u5)
    {
        if (u1 == 0)
        {
            rho = 0.0;
            u = 0.0;
            v = 0.0;
            w = 0.0;
            p = 0.0;
        }
        else
        {
            rho = u1;
            u = u2 / u1;
            v = u3 / u1;
            w = u4 / u1;
            p = (G - 1) * (u5 - 0.5 * ((u2 * u2 + u3 * u3 + u4 * u4) / u1));

            if (p < 0)
            {
                std::cout << "p below zero, " << p << ", " << u << ", " << rho
                          << ", " << u1 << ", " << u2 << ", " << u3 << ", " << u4 << std::endl;
            }
            if (std::isnan(p))
            {
                std::cout << "p error, " << p << ", " << u << ", " << rho
                          << ", " << u1 << ", " << u2 << ", " << u3 << ", " << u4 << std::endl;
                return false;
            }
        }
        fixVacuum();
        return true;
    }

} // namespace fluid