#include "domain/Domain.hpp"
#include <cmath>
#include "constants.hpp"

namespace fluid
{

    void Domain::setup(const int i, const double x, const double y, const double z,
                       const double density, const State &initial)
    {
        id = i;
        nx = x * density + 2;
        ny = y * density + 2;
        nz = z * density + 2;
        nxFaces = nx - 1;
        nyFaces = ny - 1;
        nzFaces = nz - 1;
        boxDims = 1 / density;

        rho_.resize(nx * ny * nz, initial.rho);
        u_.resize(nx * ny * nz, initial.u);
        v_.resize(nx * ny * nz, initial.v);
        w_.resize(nx * ny * nz, initial.w);
        p_.resize(nx * ny * nz, initial.p);
        a_.resize(nx * ny * nz, initial.a);

        xFaces.resize(nxFaces * ny * nz);
        yFaces.resize(nx * nyFaces * nz);
        zFaces.resize(nx * ny * nzFaces);

        setGhostCellMasks();
    }

    int Domain::getGlobalIndex(const int &x, const int &y, const int &z)
    {
        return z + nz * (y + ny * x);
    }

    double &Domain::rho(const int &x, const int &y, const int &z)
    {
        return rho_[getGlobalIndex(x, y, z)];
    }

    double &Domain::u(const int &x, const int &y, const int &z)
    {
        return u_[getGlobalIndex(x, y, z)];
    }

    double &Domain::v(const int &x, const int &y, const int &z)
    {
        return v_[getGlobalIndex(x, y, z)];
    }

    double &Domain::w(const int &x, const int &y, const int &z)
    {
        return w_[getGlobalIndex(x, y, z)];
    }

    double &Domain::p(const int &x, const int &y, const int &z)
    {
        return p_[getGlobalIndex(x, y, z)];
    }

    double &Domain::a(const int &x, const int &y, const int &z)
    {
        /*if (rho(x, y, z) == 0.0)
        {
            return a_[getGlobalIndex(x, y, z)] = 0.0;
        }
        else
        {
            return a_[getGlobalIndex(x, y, z)] = sqrt(G * p(x, y, z) / rho(x, y, z));
        }*/
        return a_[getGlobalIndex(x, y, z)];
    }

    Flux &Domain::xfAt(const int &x, const int &y, const int &z)
    {
        return xFaces[x + nxFaces * (y + ny * z)];
    }

    Flux &Domain::yfAt(const int &x, const int &y, const int &z)
    {
        return yFaces[x + nx * (y + nyFaces * z)];
    }

    Flux &Domain::zfAt(const int &x, const int &y, const int &z)
    {
        return zFaces[x + nx * (y + ny * z)];
    }

    StateRef Domain::at(const int &x, const int &y, const int &z)
    {
        return StateRef{rho(x, y, z), u(x, y, z), v(x, y, z), w(x, y, z), p(x, y, z), a(x, y, z)};
    }

    void Domain::setGhostCellMasks()
    {
        ghostCellMask.resize(nx * ny * nz, 0);

        // Set masks for x-y faces
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                ghostCellMask[0 + nz * (j + ny * i)] = 1;
                ghostCellMask[(nz - 1) + nz * (j + ny * i)] = 1;
            }
        }

        // Set masks for x-z faces
        for (int i = 0; i < nx; i++)
        {
            for (int k = 0; k < nz; k++)
            {
                ghostCellMask[k + nz * (0 + ny * i)] = 1;
                ghostCellMask[k + nz * ((ny - 1) + ny * i)] = 1;
            }
        }

        // Set masks for y-z faces
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                ghostCellMask[k + nz * (j + ny * 0)] = 1;
                ghostCellMask[k + nz * (j + ny * (nx - 1))] = 1;
            }
        }
    }

} // namespace fluid