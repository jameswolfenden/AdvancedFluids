#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include <cstdint>
#include "types/State.hpp"
#include "types/Flux.hpp"

namespace fluid
{

    class Domain
    {
    public:
        int nx, ny, nz, id;
        int nxFaces, nyFaces, nzFaces;
        double xOrigin, yOrigin, zOrigin;
        double boxDims;
        Domain *sides[6];
        std::vector<double> rho_;
        std::vector<double> u_;
        std::vector<double> v_;
        std::vector<double> w_;
        std::vector<double> p_;
        std::vector<double> a_;
        std::vector<uint8_t> ghostCellMask;

        void setup(const int i, const double x, const double y, const double z,
                   const double density, const State &initial);

        int getGlobalIndex(const int &x, const int &y, const int &z);

        // Access methods for state variables
        double &rho(const int &x, const int &y, const int &z);
        double &u(const int &x, const int &y, const int &z);
        double &v(const int &x, const int &y, const int &z);
        double &w(const int &x, const int &y, const int &z);
        double &p(const int &x, const int &y, const int &z);
        double &a(const int &x, const int &y, const int &z);

        // Access methods for faces
        Flux &xfAt(const int &x, const int &y, const int &z);
        Flux &yfAt(const int &x, const int &y, const int &z);
        Flux &zfAt(const int &x, const int &y, const int &z);

        StateRef at(const int &x, const int &y, const int &z);

    private:
        std::vector<Flux> xFaces, yFaces, zFaces;

        void setGhostCellMasks();
    };

} // namespace fluid

#endif // DOMAIN_HPP