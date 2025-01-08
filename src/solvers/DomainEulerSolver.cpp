#include "solvers/DomainEulerSolver.hpp"
#include "types/ConservedState.hpp"
#include <algorithm>

namespace fluid
{

    DomainEulerSolver::DomainEulerSolver(double clf) : clf(clf), minT(0.0) {}

    bool DomainEulerSolver::updateDomains(std::vector<Domain> &domains)
    {
        if (!updateGhostCells(domains))
            return false;
        if (!fetchTimeStep(domains))
            return false;

        // First order splitting scheme
        if (!updateX(domains))
            return false;
        if (!updateGhostCells(domains))
            return false;
        if (!updateY(domains))
            return false;
        if (!updateGhostCells(domains))
            return false;
        if (!updateZ(domains))
            return false;
        if (!updateGhostCells(domains))
            return false;

        // Second order correction
        if (!updateZ(domains))
            return false;
        if (!updateGhostCells(domains))
            return false;
        if (!updateY(domains))
            return false;
        if (!updateGhostCells(domains))
            return false;
        if (!updateX(domains))
            return false;
        if (!updateGhostCells(domains))
            return false;

        return true;
    }

    bool DomainEulerSolver::fetchTimeStep(std::vector<Domain> &domains)
    {
        minT = 1e10;
        for (auto &domain : domains)
        {
            if (!timeStep(domain))
                return false;
        }
        return true;
    }

    bool DomainEulerSolver::timeStep(Domain &domain)
    {
        for (int i = 1; i < domain.nx - 1; i++)
        {
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    double local_min_t;

                    local_min_t = clf * domain.boxDims /
                                  (domain.u(i, j, k) + domain.a(i, j, k));
                    minT = std::min(minT, local_min_t);

                    local_min_t = clf * domain.boxDims /
                                  (domain.v(i, j, k) + domain.a(i, j, k));
                    minT = std::min(minT, local_min_t);

                    local_min_t = clf * domain.boxDims /
                                  (domain.w(i, j, k) + domain.a(i, j, k));
                    minT = std::min(minT, local_min_t);
                }
            }
        }
        return true;
    }

    bool DomainEulerSolver::updateX(std::vector<Domain> &domains)
    {
        for (auto &domain : domains)
        {
            if (!xFaces(domain))
                return false;
            if (!xBoxes(domain))
                return false;
        }
        return true;
    }

    bool DomainEulerSolver::updateY(std::vector<Domain> &domains)
    {
        for (auto &domain : domains)
        {
            if (!yFaces(domain))
                return false;
            if (!yBoxes(domain))
                return false;
        }
        return true;
    }

    bool DomainEulerSolver::updateZ(std::vector<Domain> &domains)
    {
        for (auto &domain : domains)
        {
            if (!zFaces(domain))
                return false;
            if (!zBoxes(domain))
                return false;
        }
        return true;
    }

    bool DomainEulerSolver::xFaces(Domain &domain)
    {
        std::atomic<bool> errorFlag(false);

#pragma omp parallel for collapse(3)
        for (int i = 0; i < domain.nxFaces; i++)
        {
            for (int j = 0; j < domain.ny; j++)
            {
                for (int k = 0; k < domain.nz; k++)
                {
                    if (errorFlag.load())
                        continue;

                    StateView left(domain.rho(i, j, k), domain.u(i, j, k),
                                   domain.v(i, j, k), domain.w(i, j, k),
                                   domain.p(i, j, k), domain.a(i, j, k));
                    StateView right(domain.rho(i + 1, j, k), domain.u(i + 1, j, k),
                                    domain.v(i + 1, j, k), domain.w(i + 1, j, k),
                                    domain.p(i + 1, j, k), domain.a(i + 1, j, k));

                    if (!rs.findStar(left, right, domain.xfAt(i, j, k)))
                    {

                        errorFlag.store(true);
                    }
                }
            }
        }
        return !errorFlag.load();
    }

    bool DomainEulerSolver::yFaces(Domain &domain)
    {
        std::atomic<bool> errorFlag(false);

#pragma omp parallel for collapse(3)
        for (int i = 0; i < domain.nx; i++)
        {
            for (int j = 0; j < domain.nyFaces; j++)
            {
                for (int k = 0; k < domain.nz; k++)
                {
                    if (errorFlag.load())
                        continue;

                    StateView left(domain.rho(i, j, k), domain.v(i, j, k),
                                   domain.u(i, j, k), domain.w(i, j, k),
                                   domain.p(i, j, k), domain.a(i, j, k));
                    StateView right(domain.rho(i, j + 1, k), domain.v(i, j + 1, k),
                                    domain.u(i, j + 1, k), domain.w(i, j + 1, k),
                                    domain.p(i, j + 1, k), domain.a(i, j + 1, k));

                    if (!rs.findStar(left, right, domain.yfAt(i, j, k)))
                    {

                        errorFlag.store(true);
                    }
                }
            }
        }
        return !errorFlag.load();
    }

    bool DomainEulerSolver::zFaces(Domain &domain)
    {
        std::atomic<bool> errorFlag(false);

#pragma omp parallel for collapse(3)
        for (int i = 0; i < domain.nx; i++)
        {
            for (int j = 0; j < domain.ny; j++)
            {
                for (int k = 0; k < domain.nzFaces; k++)
                {
                    if (errorFlag.load())
                        continue;

                    StateView left(domain.rho(i, j, k), domain.w(i, j, k),
                                   domain.v(i, j, k), domain.u(i, j, k),
                                   domain.p(i, j, k), domain.a(i, j, k));
                    StateView right(domain.rho(i, j, k + 1), domain.w(i, j, k + 1),
                                    domain.v(i, j, k + 1), domain.u(i, j, k + 1),
                                    domain.p(i, j, k + 1), domain.a(i, j, k + 1));

                    if (!rs.findStar(left, right, domain.zfAt(i, j, k)))
                    {

                        errorFlag.store(true);
                    }
                }
            }
        }
        return !errorFlag.load();
    }

    bool DomainEulerSolver::xBoxes(Domain &domain)
    {
        std::atomic<bool> errorFlag(false);

#pragma omp parallel for collapse(3)
        for (int i = 1; i < domain.nx - 1; i++)
        {
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    if (errorFlag.load())
                        continue;

                    ConservedState state(domain.rho(i, j, k), domain.u(i, j, k),
                                         domain.v(i, j, k), domain.w(i, j, k),
                                         domain.p(i, j, k), domain.a(i, j, k));

                    double u1 = state.u1() + minT * (domain.xfAt(i - 1, j, k).f1 - domain.xfAt(i, j, k).f1) / domain.boxDims;
                    double u2 = state.u2() + minT * (domain.xfAt(i - 1, j, k).f2 - domain.xfAt(i, j, k).f2) / domain.boxDims;
                    double u3 = state.u3() + minT * (domain.xfAt(i - 1, j, k).f3 - domain.xfAt(i, j, k).f3) / domain.boxDims;
                    double u4 = state.u4() + minT * (domain.xfAt(i - 1, j, k).f4 - domain.xfAt(i, j, k).f4) / domain.boxDims;
                    double u5 = state.u5() + minT * (domain.xfAt(i - 1, j, k).f5 - domain.xfAt(i, j, k).f5) / domain.boxDims;

                    if (!state.updateFromConservatives(u1, u2, u3, u4, u5))
                    {
                        errorFlag.store(true);
                    }
                }
            }
        }
        return !errorFlag.load();
    }

    bool DomainEulerSolver::yBoxes(Domain &domain)
    {
        std::atomic<bool> errorFlag(false);

#pragma omp parallel for collapse(3)
        for (int i = 1; i < domain.nx - 1; i++)
        {
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    if (errorFlag.load())
                        continue;

                    ConservedState state(domain.rho(i, j, k), domain.v(i, j, k),
                                         domain.u(i, j, k), domain.w(i, j, k),
                                         domain.p(i, j, k), domain.a(i, j, k));

                    double u1 = state.u1() + minT * (domain.yfAt(i, j - 1, k).f1 - domain.yfAt(i, j, k).f1) / domain.boxDims;
                    double u2 = state.u2() + minT * (domain.yfAt(i, j - 1, k).f2 - domain.yfAt(i, j, k).f2) / domain.boxDims;
                    double u3 = state.u3() + minT * (domain.yfAt(i, j - 1, k).f3 - domain.yfAt(i, j, k).f3) / domain.boxDims;
                    double u4 = state.u4() + minT * (domain.yfAt(i, j - 1, k).f4 - domain.yfAt(i, j, k).f4) / domain.boxDims;
                    double u5 = state.u5() + minT * (domain.yfAt(i, j - 1, k).f5 - domain.yfAt(i, j, k).f5) / domain.boxDims;

                    if (!state.updateFromConservatives(u1, u2, u3, u4, u5))
                    {
                        errorFlag.store(true);
                    }
                }
            }
        }
        return !errorFlag.load();
    }

    bool DomainEulerSolver::zBoxes(Domain &domain)
    {
        std::atomic<bool> errorFlag(false);

#pragma omp parallel for collapse(3)
        for (int i = 1; i < domain.nx - 1; i++)
        {
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    if (errorFlag.load())
                        continue;

                    ConservedState state(domain.rho(i, j, k), domain.w(i, j, k),
                                         domain.v(i, j, k), domain.u(i, j, k),
                                         domain.p(i, j, k), domain.a(i, j, k));

                    double u1 = state.u1() + minT * (domain.zfAt(i, j, k - 1).f1 - domain.zfAt(i, j, k).f1) / domain.boxDims;
                    double u2 = state.u2() + minT * (domain.zfAt(i, j, k - 1).f2 - domain.zfAt(i, j, k).f2) / domain.boxDims;
                    double u3 = state.u3() + minT * (domain.zfAt(i, j, k - 1).f3 - domain.zfAt(i, j, k).f3) / domain.boxDims;
                    double u4 = state.u4() + minT * (domain.zfAt(i, j, k - 1).f4 - domain.zfAt(i, j, k).f4) / domain.boxDims;
                    double u5 = state.u5() + minT * (domain.zfAt(i, j, k - 1).f5 - domain.zfAt(i, j, k).f5) / domain.boxDims;

                    if (!state.updateFromConservatives(u1, u2, u3, u4, u5))
                    {
                        errorFlag.store(true);
                    }
                }
            }
        }
        return !errorFlag.load();
    }

    bool DomainEulerSolver::updateGhostCells(std::vector<Domain> &domains)
    {
        for (auto &domain : domains)
        {
            if (!xGhosts(domain))
                return false;
            if (!yGhosts(domain))
                return false;
            if (!zGhosts(domain))
                return false;
        }
        return true;
    }

    bool DomainEulerSolver::xGhosts(Domain &domain)
    {
        // Handle +x side (side 0)
        if (domain.sides[0])
        {
            // Transmissive boundary condition
#pragma omp parallel for collapse(2)
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    domain.at(domain.nx - 1, j, k) = domain.sides[0]->at(1, j, k);
                }
            }
        }
        else
        {
            // Reflective boundary condition
#pragma omp parallel for collapse(2)
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    domain.at(domain.nx - 1, j, k) = domain.at(domain.nx - 2, j, k);
                    domain.u(domain.nx - 1, j, k) = -domain.u(domain.nx - 1, j, k);
                }
            }
        }

        // Handle -x side (side 1)
        if (domain.sides[1])
        {
            // Transmissive boundary condition
#pragma omp parallel for collapse(2)
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    domain.at(0, j, k) = domain.sides[1]->at(domain.sides[1]->nx - 2, j, k);
                }
            }
        }
        else
        {
            // Reflective boundary condition
#pragma omp parallel for collapse(2)
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    domain.at(0, j, k) = domain.at(1, j, k);
                    domain.u(0, j, k) = -domain.u(0, j, k);
                }
            }
        }
        return true;
    }

    bool DomainEulerSolver::yGhosts(Domain &domain)
    {
        // Handle +y side (side 2)
        if (domain.sides[2])
        {
            // Transmissive boundary condition
#pragma omp parallel for collapse(2)
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    domain.at(i, domain.ny - 1, k) = domain.sides[2]->at(i, 1, k);
                }
            }
        }
        else
        {
            // Reflective boundary condition
#pragma omp parallel for collapse(2)
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    domain.at(i, domain.ny - 1, k) = domain.at(i, domain.ny - 2, k);
                    domain.v(i, domain.ny - 1, k) = -domain.v(i, domain.ny - 1, k);
                }
            }
        }

        // Handle -y side (side 3)
        if (domain.sides[3])
        {
            // Transmissive boundary condition
#pragma omp parallel for collapse(2)
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    domain.at(i, 0, k) = domain.sides[3]->at(i, domain.sides[3]->ny - 2, k);
                }
            }
        }
        else
        {
            // Reflective boundary condition
#pragma omp parallel for collapse(2)
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    domain.at(i, 0, k) = domain.at(i, 1, k);
                    domain.v(i, 0, k) = -domain.v(i, 0, k);
                }
            }
        }
        return true;
    }

    bool DomainEulerSolver::zGhosts(Domain &domain)
    {
        // Handle +z side (side 4)
        if (domain.sides[4])
        {
            // Transmissive boundary condition
#pragma omp parallel for collapse(2)
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int j = 1; j < domain.ny - 1; j++)
                {
                    domain.at(i, j, domain.nz - 1) = domain.sides[4]->at(i, j, 1);
                }
            }
        }
        else
        {
            // Reflective boundary condition
#pragma omp parallel for collapse(2)
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int j = 1; j < domain.ny - 1; j++)
                {
                    domain.at(i, j, domain.nz - 1) = domain.at(i, j, domain.nz - 2);
                    domain.w(i, j, domain.nz - 1) = -domain.w(i, j, domain.nz - 1);
                }
            }
        }

        // Handle -z side (side 5)
        if (domain.sides[5])
        {
            // Transmissive boundary condition
#pragma omp parallel for collapse(2)
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int j = 1; j < domain.ny - 1; j++)
                {
                    domain.at(i, j, 0) = domain.sides[5]->at(i, j, domain.sides[5]->nz - 2);
                }
            }
        }
        else
        {
            // Reflective boundary condition
#pragma omp parallel for collapse(2)
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int j = 1; j < domain.ny - 1; j++)
                {
                    domain.at(i, j, 0) = domain.at(i, j, 1);
                    domain.w(i, j, 0) = -domain.w(i, j, 0);
                }
            }
        }
        return true;
    }

} // namespace fluid