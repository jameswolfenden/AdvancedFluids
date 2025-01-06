#ifndef RIEMANN_SOLVER_HPP
#define RIEMANN_SOLVER_HPP

#include "../types/Flux.hpp"
#include "../constants.hpp"

namespace fluid
{

    class RiemannSolver
    {
    public:
        bool findStar(const double &rhoL, const double &uL, const double &vL, const double &wL,
                      const double &aL, const double &pL, const double &rhoR, const double &uR,
                      const double &vR, const double &wR, const double &aR, const double &pR,
                      Flux &fl);

    private:
        static constexpr double TOL = 0.000001;
        static constexpr double G1 = (G - 1) / (2 * G);
        static constexpr double G2 = (G + 1) / (2 * G);
        static constexpr double G3 = 2 * G / (G - 1);
        static constexpr double G4 = 2 / (G - 1);
        static constexpr double G5 = 2 / (G + 1);
        static constexpr double G6 = (G - 1) / (G + 1);
        static constexpr double G7 = (G - 1) / 2;
        static constexpr double G8 = G - 1;

        bool pickSide(const double &rhoL, const double &vL, const double &wL, const double &pL,
                      const double &rhoR, const double &vR, const double &wR, const double &pR,
                      Flux &fl, double &tempP, double &tempU);

        bool testVacuum(const double &rhoL, const double &uL, const double &vL, const double &wL,
                        const double &aL, const double &pL, const double &rhoR, const double &uR,
                        const double &vR, const double &wR, const double &aR, const double &pR,
                        Flux &fl);

        bool pickStartVal(const int errorStage, const double &rhoL, const double &uL,
                          const double &aL, const double &pL, const double &rhoR, const double &uR,
                          const double &aR, const double &pR, double &tempP);

        bool iterateP(const double &rhoL, const double &uL, const double &aL, const double &pL,
                      const double &rhoR, const double &uR, const double &aR, const double &pR,
                      double &tempP, double &tempU);
    };

} // namespace fluid

#endif // RIEMANN_SOLVER_HPP