#include "types/State.hpp"

namespace fluid
{
    State::State(double rho, double u, double v, double w, double p)
        : rho(rho), u(u), v(v), w(w), p(p)
    {
        a = sqrt(G * p / rho); // Calculate speed of sound
    }

    // Copy constructor, allows for copying of state objects using the '=' operator
    StateRef &StateRef::operator=(const StateRef &other)
    {
        if (this != &other)
        {
            // If the reference is to a different state, update the state from the other state
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