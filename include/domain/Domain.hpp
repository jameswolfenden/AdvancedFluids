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

        int getGlobalIndex(int x, int y, int z);

        void setup(const int i, const double x, const double y, const double z,
                   const double density, const State &initial);

        // Access methods for state variables
        double &rho(int x, int y, int z);
        double &u(int x, int y, int z);
        double &v(int x, int y, int z);
        double &w(int x, int y, int z);
        double &p(int x, int y, int z);
        double &a(int x, int y, int z);

        // Access methods for faces
        Flux &xfAt(int x, int y, int z);
        Flux &yfAt(int x, int y, int z);
        Flux &zfAt(int x, int y, int z);

        StateRef at(int x, int y, int z);

    private:
        std::vector<Flux> xFaces, yFaces, zFaces;

        void setGhostCellMasks();
    };

} // namespace fluid

#endif // DOMAIN_HPP