#ifndef STATE_HPP
#define STATE_HPP

#include <cmath>
#include "constants.hpp"

namespace fluid
{
    // Struct to hold the state of the fluid, used mostly to initialise the properties of a domain
    struct State
    {
        double rho, u, v, w, p, a;

        State(double rho, double u, double v, double w, double p);
    };

    // Struct to hold references to the properties a box in the domain to allow for easy access and modification
    struct StateRef
    {
        double &rho, &u, &v, &w, &p, &a;

        StateRef(double &rho, double &u, double &v, double &w, double &p, double &a) : rho(rho), u(u), v(v), w(w), p(p), a(a) {}
        explicit StateRef(State &state) : rho(state.rho), u(state.u), v(state.v), w(state.w), p(state.p), a(state.a) {}

        StateRef &operator=(const StateRef &other);

    };

} // namespace fluid

#endif // STATE_HPP