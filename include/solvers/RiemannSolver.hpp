#ifndef RIEMANN_SOLVER_HPP
#define RIEMANN_SOLVER_HPP

#include "../types/Flux.hpp"
#include "../constants.hpp"
#include "../types/State.hpp"

namespace fluid
{

    class RiemannSolver
    {
    public:
        explicit RiemannSolver(std::pair<double, int> vals = {1E-6, 10000}) : TOL(vals.first), MAX_ITER(vals.second) {}

        bool findStar(const StateView &left, const StateView &right, Flux &fl);

        double getLastP() const { return lastP; }

    private:
        static constexpr double G1 = (G - 1) / (2 * G);
        static constexpr double G2 = (G + 1) / (2 * G);
        static constexpr double G3 = 2 * G / (G - 1);
        static constexpr double G4 = 2 / (G - 1);
        static constexpr double G5 = 2 / (G + 1);
        static constexpr double G6 = (G - 1) / (G + 1);
        static constexpr double G7 = (G - 1) / 2;
        static constexpr double G8 = G - 1;

        double lastP = -1.0;
        double TOL;
        int MAX_ITER;

        bool pickSide(const StateView &left, const StateView &right,
                      Flux &fl, double &tempP, double &tempU);

        bool testVacuum(const StateView &left, const StateView &right, Flux &fl);

        bool pickStartVal(const int errorStage, const StateView &left, const StateView &right, double &tempP);

        bool iterateP(const StateView &left, const StateView &right,
                      double &tempP, double &tempU);

        class Impl;
    };

} // namespace fluid

#endif // RIEMANN_SOLVER_HPP