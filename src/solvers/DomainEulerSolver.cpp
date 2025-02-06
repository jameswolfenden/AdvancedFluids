#include "solvers/DomainEulerSolver.hpp"
#include "types/ConservedState.hpp"

namespace fluid
{

    DomainEulerSolver::DomainEulerSolver(double clf) : clf(clf), minT(0.0) {}

    bool DomainEulerSolver::updateDomains(std::vector<Domain> &domains)
    {
        // Update ghost cells first for safety
        updateGhostCells(domains);
        fetchTimeStep(domains);

        // First order splitting scheme (see Toro, 2009)
        if (!updateX(domains))
            return false;
        updateGhostCells(domains);
        if (!updateY(domains))
            return false;
        updateGhostCells(domains);
        if (!updateZ(domains))
            return false;
        updateGhostCells(domains);

        // Second order correction
        if (!updateZ(domains))
            return false;
        updateGhostCells(domains);
        if (!updateY(domains))
            return false;
        updateGhostCells(domains);
        if (!updateX(domains))
            return false;
        updateGhostCells(domains);

        return true;
    }

    void DomainEulerSolver::fetchTimeStep(std::vector<Domain> &domains)
    {
        minT = 1e10;
        for (auto &domain : domains)
        {
            timeStep(domain);
        }
    }

    void DomainEulerSolver::timeStep(Domain &domain)
    {
        // Only runs once per timestep but could be optimised
        // Checks the CFL condition for each cell and returns the minimum
        for (int i = 1; i < domain.nx - 1; i++)
        {
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    double local_min_t = clf * domain.boxDims /
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
    }

    bool DomainEulerSolver::updateX(std::vector<Domain> &domains)
    {
        std::atomic<bool> errorFlag(false); // Atomic for thread safety
        #pragma omp parallel for // Parallelise the loop, see the parallelDomainEulerSolver benchmark for performance comparison
        for (auto &domain : domains)
        {
            if (!xFaces(domain))
                errorFlag.store(true);
            if (!xBoxes(domain))
                errorFlag.store(true);
        }
        return !errorFlag.load();
    }

    bool DomainEulerSolver::updateY(std::vector<Domain> &domains)
    {
        std::atomic<bool> errorFlag(false); // Atomic for thread safety
        #pragma omp parallel for // Parallelise the loop
        for (auto &domain : domains)
        {
            if (!yFaces(domain))
                errorFlag.store(true);
            if (!yBoxes(domain))
                errorFlag.store(true);
        }
        return !errorFlag.load();
    }

    bool DomainEulerSolver::updateZ(std::vector<Domain> &domains)
    {
        std::atomic<bool> errorFlag(false); // Atomic for thread safety
        #pragma omp parallel for // Parallelise the loop
        for (auto &domain : domains)
        {
            if (!zFaces(domain))
                errorFlag.store(true);
            if (!zBoxes(domain))
                errorFlag.store(true);
        }
        return !errorFlag.load();
    }

    bool DomainEulerSolver::xFaces(Domain &domain)
    {
        std::atomic<bool> errorFlag(false);
        for (int i = 0; i < domain.nxFaces; i++)
        {
            for (int j = 0; j < domain.ny; j++)
            {
                for (int k = 0; k < domain.nz; k++)
                {
                    if (errorFlag.load())
                        continue; // Keep skipping if an error has been found

                    int index = k + domain.nz * (j + domain.ny * i);

                    StateRef left(domain.rho_[index], domain.u_[index], domain.v_[index],
                                   domain.w_[index], domain.p_[index], domain.a_[index]);

                    index += domain.nz * domain.ny;

                    StateRef right(domain.rho_[index], domain.u_[index], domain.v_[index],
                                    domain.w_[index], domain.p_[index], domain.a_[index]);

                    // Find the state at the x face from the left and right states
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
        for (int i = 0; i < domain.nx; i++)
        {
            for (int j = 0; j < domain.nyFaces; j++)
            {
                for (int k = 0; k < domain.nz; k++)
                {
                    if (errorFlag.load())
                        continue; // Keep skipping if an error has been found

                    int index = k + domain.nz * (j + domain.ny * i);

                    // Note the order of the arguments is different to xFaces so the velocities are transformed onto the correct axis for 'left' and 'right'
                    StateRef left(domain.rho_[index], domain.v_[index], domain.u_[index],
                                   domain.w_[index], domain.p_[index], domain.a_[index]);

                    index += domain.nz;

                    StateRef right(domain.rho_[index], domain.v_[index], domain.u_[index],
                                    domain.w_[index], domain.p_[index], domain.a_[index]);

                    // Find the state at the y face from the left and right states
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
        for (int i = 0; i < domain.nx; i++)
        {
            for (int j = 0; j < domain.ny; j++)
            {
                for (int k = 0; k < domain.nzFaces; k++)
                {
                    if (errorFlag.load())
                        continue; // Keep skipping if an error has been found

                    int index = k + domain.nz * (j + domain.ny * i);

                    // Note the order of the arguments is different to xFaces so the velocities are transformed onto the correct axis for 'left' and 'right'
                    StateRef left(domain.rho_[index], domain.w_[index], domain.v_[index],
                                   domain.u_[index], domain.p_[index], domain.a_[index]);

                    index++;

                    StateRef right(domain.rho_[index], domain.w_[index], domain.v_[index],
                                    domain.u_[index], domain.p_[index], domain.a_[index]);

                    // Find the state at the z face from the left and right states
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
        for (int i = 1; i < domain.nx - 1; i++)
        {
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    if (errorFlag.load())
                        continue;

                    // Conserved State is a helper class to convert between primitive and conserved variables
                    // Conserved variables are needed for updating based on fluxes generated by the Riemann solver
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
        for (int i = 1; i < domain.nx - 1; i++)
        {
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    if (errorFlag.load())
                        continue;

                    // Note that again the order of the arguments is different to xBoxes to transform the velocities back to the correct axis
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
        for (int i = 1; i < domain.nx - 1; i++)
        {
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    if (errorFlag.load())
                        continue;

                    // Note that again the order of the arguments is different to xBoxes to transform the velocities back to the correct axis
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

    void DomainEulerSolver::updateGhostCells(std::vector<Domain> &domains)
    {
        #pragma omp parallel for // Do each domain in parallel
        for (auto &domain : domains)
        {
            xGhosts(domain);
            yGhosts(domain);
            zGhosts(domain);
        }
    }

    void DomainEulerSolver::xGhosts(Domain &domain)
    {
        // Handle +x side (side 0)
        if (domain.sides[0])
        {
            // Transmissive boundary condition
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
            for (int j = 1; j < domain.ny - 1; j++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    domain.at(0, j, k) = domain.at(1, j, k);
                    domain.u(0, j, k) = -domain.u(0, j, k);
                }
            }
        }
    }

    void DomainEulerSolver::yGhosts(Domain &domain)
    {
        // Handle +y side (side 2)
        if (domain.sides[2])
        {
            // Transmissive boundary condition
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
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int k = 1; k < domain.nz - 1; k++)
                {
                    domain.at(i, 0, k) = domain.at(i, 1, k);
                    domain.v(i, 0, k) = -domain.v(i, 0, k);
                }
            }
        }
    }

    void DomainEulerSolver::zGhosts(Domain &domain)
    {
        // Handle +z side (side 4)
        if (domain.sides[4])
        {
            // Transmissive boundary condition
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
            for (int i = 1; i < domain.nx - 1; i++)
            {
                for (int j = 1; j < domain.ny - 1; j++)
                {
                    domain.at(i, j, 0) = domain.at(i, j, 1);
                    domain.w(i, j, 0) = -domain.w(i, j, 0);
                }
            }
        }
    }

} // namespace fluid