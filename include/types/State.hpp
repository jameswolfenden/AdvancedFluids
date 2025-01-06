#ifndef STATE_HPP
#define STATE_HPP

#include <cmath>
#include "constants.hpp"

namespace fluid
{

    struct State
    {
        double rho, u, v, w, p, a;

        State();
        State(double rho, double u, double v, double w, double p);
    };

    struct StateRef
    {
        double &rho, &u, &v, &w, &p, &a;

        StateRef &operator=(const StateRef &other);
    };

} // namespace fluid

#endif // STATE_HPP