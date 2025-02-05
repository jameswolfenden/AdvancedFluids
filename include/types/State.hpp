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

        State();
        State(double rho, double u, double v, double w, double p);
    };

    // Struct to hold references to the properties a box in the domain to allow for easy access and modification
    struct StateView
    {
        double &rho, &u, &v, &w, &p, &a;

        StateView(double &rho, double &u, double &v, double &w, double &p, double &a) : rho(rho), u(u), v(v), w(w), p(p), a(a) {}
        explicit StateView(State &state) : rho(state.rho), u(state.u), v(state.v), w(state.w), p(state.p), a(state.a) {}
    };
    
    // Allows for the properties of one box in the domain to be copied to another
    struct StateRef
    {
        double &rho, &u, &v, &w, &p, &a;

        StateRef &operator=(const StateRef &other);
    };

} // namespace fluid

#endif // STATE_HPP