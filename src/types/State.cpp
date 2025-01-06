#include "types/State.hpp"

namespace fluid
{

    State::State() = default;

    State::State(double rho, double u, double v, double w, double p)
        : rho(rho), u(u), v(v), w(w), p(p)
    {
        a = sqrt(G * p / rho);
    }

    StateRef &StateRef::operator=(const StateRef &other)
    {
        if (this != &other)
        {
            rho = other.rho;
            u = other.u;
            v = other.v;
            w = other.w;
            p = other.p;
            a = other.a;
        }
        return *this;
    }

} // namespace fluid