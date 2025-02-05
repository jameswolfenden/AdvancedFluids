#ifndef DOMAIN_EULER_SOLVER_HPP
#define DOMAIN_EULER_SOLVER_HPP

#include <vector>
#include <atomic>
#include "domain/Domain.hpp"
#include "RiemannSolver.hpp"

namespace fluid
{

    class DomainEulerSolverBenchmark;

    class DomainEulerSolver
    {
    public:
        explicit DomainEulerSolver(double clf);

        bool updateDomains(std::vector<Domain> &domains);
        double getMinTimeStep() const { return minT; }

    private:
        friend class DomainEulerSolverBenchmark; // Allows the benchmark to access private members
        RiemannSolver rs;
        double clf;
        double minT;

        // Time step calculation
        void fetchTimeStep(std::vector<Domain> &domains);
        void timeStep(Domain &domain);

        // Domain updates
        bool updateX(std::vector<Domain> &domains);
        bool updateY(std::vector<Domain> &domains);
        bool updateZ(std::vector<Domain> &domains);

        // Face calculations
        bool xFaces(Domain &domain);
        bool yFaces(Domain &domain);
        bool zFaces(Domain &domain);

        // Box updates
        bool xBoxes(Domain &domain);
        bool yBoxes(Domain &domain);
        bool zBoxes(Domain &domain);

        // Ghost cell handling
        void updateGhostCells(std::vector<Domain> &domains);
        void xGhosts(Domain &domain);
        void yGhosts(Domain &domain);
        void zGhosts(Domain &domain);
    };

} // namespace fluid

#endif // DOMAIN_EULER_SOLVER_HPP