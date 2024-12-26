#include <cmath>
#include <iostream>
#include <vector>
#include "H5Cpp.h"
#include <fstream>

const double G = 1.4;

struct State
{
    double rho, u, v, w, p, a;
    State() {}
    State(double rho, double u, double v, double w, double p) : rho(rho), u(u), v(v), w(w), p(p) { a = sqrt(G * p / rho); }
};

struct StateRef
{
    double &rho, &u, &v, &w, &p, &a;
    StateRef &operator=(const StateRef &other)
    {
        if (this != &other)
        {
            rho = other.rho;
            u = other.u;
            v = other.v;
            w = other.w;
            p = other.p;
            a = other.a;
        }
        return *this;
    }
};

struct Flux
{
    double f1, f2, f3, f4, f5;
    Flux() : f1(0.0), f2(0.0), f3(0.0), f4(0.0), f5(0.0) {}
    Flux(double f1, double f2, double f3, double f4, double f5) : f1(f1), f2(f2), f3(f3), f4(f4), f5(f5) {}
    void updateFromPrimatives(double rho, double u, double v, double w, double p)
    {
        f1 = rho * u;
        f2 = rho * u * u + p;
        f3 = rho * v * u;
        f4 = rho * w * u;
        f5 = u * (rho * (0.5 * (u * u + v * v + w * w) + p / ((G - 1) * rho) + p));
    }
    // Flux(double rho, double u, double v, double w, double p) : f1(rho*u), f2(rho*u*u + p), f3(rho*v*u), f4(rho*w*u), f5(u*(0.5*(u*u + v*v + w*w) + p/((G - 1)*rho) + p)) {}
};

class RiemannSolver
{
    // check if these should be static
    const double TOL = 0.000001;
    const double G1 = (G - 1) / (2 * G);
    const double G2 = (G + 1) / (2 * G);
    const double G3 = 2 * G / (G - 1);
    const double G4 = 2 / (G - 1);
    const double G5 = 2 / (G + 1);
    const double G6 = (G - 1) / (G + 1);
    const double G7 = (G - 1) / 2;
    const double G8 = G - 1;

    bool pickSide(const double &rhoL, const double &vL, const double &wL, const double &pL, const double &rhoR, const double &vR, const double &wR, const double &pR, Flux &fl, double &tempP, double &tempU)
    {
        // pick rho value depending on the side of the discontinuity
        if (tempU >= 0.0) // pick left side
        {
            if (tempP > pL)
                fl.updateFromPrimatives(rhoL * (((tempP / pL) + G6) / (G6 * (tempP / pL) + 1)), tempU, vL, wL, tempP);
            else
                fl.updateFromPrimatives(rhoL * pow((tempP / pL), 1 / G), tempU, vL, wL, tempP);
        }
        else // pick right side
        {
            if (tempP > pR)
                fl.updateFromPrimatives(rhoR * (((tempP / pR) + G6) / (G6 * (tempP / pR) + 1)), tempU, vR, wR, tempP);
            else
                fl.updateFromPrimatives(rhoR * pow((tempP / pR), 1 / G), tempU, vR, wR, tempP);
        }
        // if (std::isnan(f.rho))
        //{
        //     std::cout << "rhoerror u: " << f.u << ", rhoL: " << rhoL << ", rhoR: " << rhoR << ", vL: " << vL << ", vR: " << vR << ", pL: " << pL << ", pR: " << pR << ", p:" << f.p << std::endl;
        //     return false;
        // }
        // calca();
        return true;
    }
    bool testVacuum(const double &rhoL, const double &uL, const double &vL, const double &wL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &wR, const double &aR, const double &pR, Flux &fl)
    {
        if ((rhoL == 0.0 || pL == 0.0) || (rhoR == 0.0 || pR == 0.0))
        {
            if (((rhoL == 0.0 || pL == 0.0) && (rhoR != 0.0 || pR != 0.0)) || (((rhoR != 0.0 || pR != 0.0) && (rhoL != 0.0 || pL != 0.0)) && 0 >= (uR - aR * G4))) // vacuum left not right
            {
                if (0 >= (uR + aR))
                {
                    std::cout << "W_R" << std::endl;
                    fl.updateFromPrimatives(rhoR, uR, vR, wR, pR);
                    // calca();
                }
                else if (0 <= (uR - aR * G4))
                {
                    std::cout << "W_L" << std::endl;
                    fl.updateFromPrimatives(rhoL, uL, vL, wL, pL);
                }
                else
                {
                    std::cout << "W_RFan" << std::endl;
                    fl.updateFromPrimatives(rhoR * pow(G5 - G6 * uR / aR, G4), G5 * (-aR + G7 * uR), vR * pow(G5 - G6 * uR / aR, G4), wR * pow(G5 - G6 * uR / aR, G4), pR * pow(G5 - G6 * uR / aR, G3));
                    // calca();
                }
                return true;
            }
            else if (((rhoR == 0.0 || pR == 0.0) && (rhoL != 0.0 || pL != 0.0)) || (((rhoR != 0.0 || pR != 0.0) && (rhoL != 0.0 || pL != 0.0)) && 0 <= (uL + aL * G4))) // vacuum right not left
            {
                if (0 <= (uL - aL))
                {
                    std::cout << "W_L" << std::endl;
                    fl.updateFromPrimatives(rhoL, uL, vL, wL, pL);
                    // calca();
                }
                else if (0 >= (uL + aL * G4))
                {
                    std::cout << "W_R" << std::endl;
                    fl.updateFromPrimatives(rhoR, uR, vR, wR, pR);
                }
                else
                {
                    std::cout << "W_LFan" << std::endl;
                    fl.updateFromPrimatives(rhoL * pow(G5 + G6 * uL / aL, G4), G5 * (aL + G7 * uL), vL * pow(G5 + G6 * uL / aL, G4), wL * pow(G5 + G6 * uL / aL, G4), pL * pow(G5 + G6 * uL / aL, G3));
                    // calca();
                }
                return true;
            }
            else if (((rhoR == 0.0 || pR == 0.0) && (rhoL == 0.0 || pL == 0.0)) || ((uL + aL * G4) < 0 && (uR - aR * G4) > 0)) // vacuum left and right
            {
                std::cout << "W_0" << std::endl;
                fl.updateFromPrimatives(0.0, 0.5 * (uL + uR), 0.5 * (vL + vR), 0.5 * (wL + wR), 0.0);
                return true;
            }
        }
        return false;
    }
    bool pickStartVal(const int errorStage, const double &rhoL, const double &uL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &aR, const double &pR, double &tempP)
    {
        double p_PV = 0.5 * (pL + pR) + 0.5 * (uL - uR) * 0.25 * (rhoL + rhoR) * (aL + aR);
        p_PV = std::max(0.0, p_PV);

        // different guesses for p
        if (errorStage == 0)
        {
            double p_min = std::min(pL, pR);
            double p_max = std::max(pL, pR);
            double q_max = p_max / p_min;
            if (q_max < 2.0 && (p_min < p_PV && p_PV < p_max))
            {
                // Select PVRS Riemann solver
                tempP = p_PV;
            }
            else if (p_PV < p_min)
            {
                // Select Two-Rarefaction Riemann solver
                double p_q = pow(pL / pR, G1);
                double u_m = (p_q * uL / aL + uR / aR + G4 * (p_q - 1.0)) / (p_q / aL + 1 / aR);
                double p_TL = 1 + G7 * (uL - u_m) / aL;
                double p_TR = 1 + G7 * (u_m - uR) / aR;
                tempP = 0.5 * (pL * pow(p_TL, G3) + pR * pow(p_TR, G3));
            }
            else
            {
                // Select Two-Shock Riemann solver with
                // PVRS as estimate
                double ge_L = sqrt((G5 / rhoL) / (G6 * pL + p_PV));
                double ge_R = sqrt((G5 / rhoR) / (G6 * pR + p_PV));
                tempP = (ge_L * pL + ge_R * pR - (uR - uL)) / (ge_L + ge_R);
            }
        }
        else if (errorStage == 1)
        {
            std::cout << "went to error stage 1" << std::endl;
            tempP = pow((aL + aR - 0.5 * G8 * (uR - uL)) / (aL / pow(pL, G1) + aR / pow(pR, G1)), G3);
        }
        else if (errorStage == 2)
        {
            tempP = 0.5 * (pL + pR);
        }
        else if (errorStage == 3)
        {
            double p_PV = 0.5 * (pL + pR) + 0.5 * (uL - uR) * 0.5 * (rhoL + rhoR) * 0.5 * (aL + aR);
            // if (p_PV > TOL)
            tempP = p_PV;
            //  else
            //     p = TOL;
        }
        else if (errorStage == 4)
        {
            double p_PV = 0.5 * (pL + pR) + 0.5 * (uL - uR) * 0.5 * (rhoL + rhoR) * 0.5 * (aL + aR);
            if (p_PV > TOL)
                tempP = p_PV;
            else
                tempP = TOL;
            double AL = 2 / ((G + 1) * rhoL);
            double BL = pL * (G - 1) / (G + 1);
            double AR = 2 / ((G + 1) * rhoR);
            double BR = pR * (G - 1) / (G + 1);
            double gL = pow(AL / (tempP + BL), 0.5);
            double gR = pow(AR / (tempP + BR), 0.5);
            double p_TS = (gL * pL + gR * pR - (uR - uL)) / (gL + gR);
            //   if (p_TS > TOL)
            tempP = p_TS;
            //  else
            //      p = TOL;
        }
        else if (errorStage == 5)
        {
            tempP = p_PV;
        }
        else if (errorStage == 6)
        {
            // Select Two-Rarefaction Riemann solver
            double p_q = pow(pL / pR, G1);
            double u_m = (p_q * uL / aL + uR / aR + G4 * (p_q - 1.0)) / (p_q / aL + 1 / aR);
            double p_TL = 1 + G7 * (uL - u_m) / aL;
            double p_TR = 1 + G7 * (u_m - uR) / aR;
            tempP = 0.5 * (pL * pow(p_TL, G3) + pR * pow(p_TR, G3));
        }
        else if (errorStage == 7)
        {
            // Select Two-Shock Riemann solver with
            // PVRS as estimate
            double ge_L = sqrt((G5 / rhoL) / (G6 * pL + p_PV));
            double ge_R = sqrt((G5 / rhoR) / (G6 * pR + p_PV));
            tempP = (ge_L * pL + ge_R * pR - (uR - uL)) / (ge_L + ge_R);
        }
        else if (errorStage == 8)
        {
            tempP = 1 / (rhoL * aL + rhoR * aR) * (rhoR * aR * pL + rhoL * aL * pR + rhoL * aL * rhoR * aR * (uL - aR));
        }
        else if (errorStage == 9)
        {
            tempP = 1 * pow(10, -6);
        }
        else
        {
            std::cout << "Error converging on p in x" << std::endl;
            std::cout << "pL: " << pL << ", rhoL: " << rhoL << ", uL: " << uL << ", aL: " << aL << ", pR: " << pR << ", rhoR: " << rhoR << ", uR: " << uR << ", aR: " << aR << std::endl;
            return false; // this doesnt actually stop the program but if you see that in the console the timestep is probably too small
        }
        return true;
    }
    bool iterateP(const double &rhoL, const double &uL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &aR, const double &pR, double &tempP, double &tempU)
    {
        bool iterate = true;
        double fsL, fsR, d_fsL, d_fsR, change;
        int errorStage = 0;

        while (iterate) // loop to try all the initial p values
        {
            iterate = pickStartVal(errorStage, rhoL, uL, aL, pL, rhoR, uR, aR, pR, tempP);
            if (!iterate)
                return false;
            int count = 0;
            // std::cout << "p guessed: " << p<< std::endl;
            while (iterate) // loop to iterate the p value
            {
                double A1 = 2 / ((G + 1) * rhoL);
                double B1 = pL * (G - 1) / (G + 1);
                double A2 = 2 / ((G + 1) * rhoR);
                double B2 = pR * (G - 1) / (G + 1);
                double p1 = A1 / (tempP + B1);
                double p2 = tempP / pL;
                double p3 = p2;
                double p4 = tempP / pR;
                double p5 = p4;
                double p6 = A2 / (tempP + B2);

                if (p1 < 0 || p2 < 0 || p3 < 0 || p4 < 0 || p5 < 0 || p6 < 0)
                {
                    errorStage++;
                    break; // exit the iterate loop to try next p
                }

                p1 = pow(p1, 0.5);
                p2 = pow(p2, G1);
                p3 = pow(p3, -G2);
                p4 = pow(p4, G1);
                p5 = pow(p5, -G2);
                p6 = pow(p6, 0.5);

                if (tempP > pL) // shock
                {
                    fsL = (tempP - pL) * p1;
                    d_fsL = p1 * (1 - (tempP - pL) / (2 * (B1 + tempP)));
                }
                else // expansion
                {
                    fsL = aL * G4 * (p2 - 1);
                    d_fsL = 1 / (rhoL * aL) * p3;
                }

                if (tempP > pR) // shock
                {
                    fsR = (tempP - pR) * p6;
                    d_fsR = p6 * (1 - (tempP - pR) / (2 * (B2 + tempP)));
                }
                else // expansion
                {
                    fsR = aR * G4 * (p4 - 1);
                    d_fsR = 1 / (rhoR * aR) * p5;
                }

                double f_ = fsL + fsR - uL + uR;
                double d_f = d_fsL + d_fsR;
                change = f_ / d_f;
                //    std::cout << f << ", " << d_f << ", " << change << std::endl;
                tempP = tempP - change; // Update new estimate of p*
                count++;

                if (TOL >= 2 * fabs(change / (change + 2 * tempP))) // iteration limit (slightly different to notes as abs of entire rhs)
                {
                    iterate = false; // we have converged
                }
                if (count > 10000)
                {
                    errorStage++;
                    break;
                }
            }
        }
        tempU = 0.5 * (uL + uR) + 0.5 * (fsR - fsL); // u*
        return true;
    }

public:
    bool findStar(const double &rhoL, const double &uL, const double &vL, const double &wL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &wR, const double &aR, const double &pR, Flux &fl) // find the values at the faces between 2 cells adjacent in x
    {
        if (testVacuum(rhoL, uL, vL, wL, aL, pL, rhoR, uR, vR, wR, aR, pR, fl))
        {
            std::cout << "vacuum" << std::endl;
            return true; // a vacuum is generated, values found without iteration required NEED TO CHECK SPEED OF SOUND!!!!!
        }
        double tempP;
        double tempU;
        if (!iterateP(rhoL, uL, aL, pL, rhoR, uR, aR, pR, tempP, tempU))
        {
            std::cout << "iteration failed" << std::endl;
            return false; // iteration failed, abort
        }

        if (!pickSide(rhoL, vL, wL, pL, rhoR, vR, wR, pR, fl, tempP, tempU))
        {
            std::cout << "pick side failed, rho is nan" << std::endl;
            return false; // rho is nan
        }
        return true;
    }
};

struct BoxConserved
{

    BoxConserved(double &rho, double &u, double &v, double &w, double &p, double &a) : rho(rho), u(u), v(v), w(w), p(p), a(a) {}

    double aCalc() // find speed of sound
    {
        if (rho == 0.0)
            a = 0.0;
        else
            a = sqrt((G * p) / rho);
        return a;
    }
    void fixVacuum()
    {
        if (rho <= 0.0)
        {
            p = 0.0;
        }
        if (p <= 0.0)
        {
            rho = 0.0;
        }
        aCalc();
    }
    // the following are to get the conservatives and fluxes from the primatives
    double u1()
    {
        return rho;
    }
    double u2()
    {
        return rho * u;
    }
    double u3()
    {
        return rho * v;
    }
    double u4()
    {
        return rho * w;
    }
    double u5()
    {
        if (rho == 0.0)
            return 0.0;
        else
            return rho * (0.5 * (u * u + v * v + w * w) + p / ((G - 1) * rho));
    }
    void updateFromConservatives(double u1, double u2, double u3, double u4, double u5) // update the primatives from new conservative values
    {
        if (u1 == 0)
        {
            rho = 0.0;
            u = 0.0;
            v = 0.0;
            w = 0.0;
            p = 0.0;
        }
        else
        {
            rho = u1;
            u = u2 / u1;
            v = u3 / u1;
            w = u4 / u1;
            p = (G - 1) * (u5 - 0.5 * ((u2 * u2 + u3 * u3 + u4 * u4) / u1));
            if (p < 0)
                std::cout << "p below zero, " << p << ", " << u << ", " << rho << ", " << u1 << ", " << u2 << ", " << u3 << ", " << u4 << std::endl;
            if (std::isnan(p))
                std::cout << "p error, " << p << ", " << u << ", " << rho << ", " << u1 << ", " << u2 << ", " << u3 << ", " << u4 << std::endl;
        }
        fixVacuum();
    }

private:
    double &rho, &u, &v, &w, &p, &a;
};

class Domain
{
public:
    int nx, ny, nz;
    int nxFaces, nyFaces, nzFaces;
    std::vector<double> rho_;
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> w_;
    std::vector<double> p_;
    std::vector<double> a_;
    std::vector<Flux> xFaces, yFaces, zFaces;
    double boxDims;
    Domain *sides[6];
    std::vector<uint8_t> ghostCellMask;
    void setup(const double x, const double y, const double z, const double density, const State &initial)
    {
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
    double &rho(const int &x, const int &y, const int &z)
    {
        return rho_[x + nx * (y + ny * z)];
    }
    double &u(const int &x, const int &y, const int &z)
    {
        return u_[x + nx * (y + ny * z)];
    }
    double &v(const int &x, const int &y, const int &z)
    {
        return v_[x + nx * (y + ny * z)];
    }
    double &w(const int &x, const int &y, const int &z)
    {
        return w_[x + nx * (y + ny * z)];
    }
    double &p(const int &x, const int &y, const int &z)
    {
        return p_[x + nx * (y + ny * z)];
    }
    double &a(const int &x, const int &y, const int &z)
    {
        if (rho(x, y, z) == 0.0)
        {
            return a_[x + nx * (y + ny * z)] = 0.0;
        }
        else
        {
            return a_[x + nx * (y + ny * z)] = sqrt(G * p(x, y, z) / rho(x, y, z));
        }
    }
    Flux &xfAt(const int &x, const int &y, const int &z)
    {
        return xFaces[x + nxFaces * (y + ny * z)];
    }
    Flux &yfAt(const int &x, const int &y, const int &z)
    {
        return yFaces[x + nx * (y + nyFaces * z)];
    }
    Flux &zfAt(const int &x, const int &y, const int &z)
    {
        return zFaces[x + nx * (y + ny * z)];
    }
    StateRef at(const int &x, const int &y, const int &z)
    {
        return StateRef{rho(x, y, z), u(x, y, z), v(x, y, z), w(x, y, z), p(x, y, z), a(x, y, z)};
    }

private:
    void setGhostCellMasks()
    {
        {
            ghostCellMask.resize(nx * ny * nz, 0);

            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    ghostCellMask[i + nx * (j + ny * 0)] = 1;
                    ghostCellMask[i + nx * (j + ny * (nz - 1))] = 1;
                }
            }
            for (int i = 0; i < nx; i++)
            {
                for (int k = 0; k < nz; k++)
                {
                    ghostCellMask[i + nx * (0 + ny * k)] = 1;
                    ghostCellMask[i + nx * ((ny - 1) + ny * k)] = 1;
                }
            }
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nz; k++)
                {
                    ghostCellMask[0 + nx * (j + ny * k)] = 1;
                    ghostCellMask[(nx - 1) + nx * (j + ny * k)] = 1;
                }
            }
        }
    }
};

class DomainSolver
{
public:
    DomainSolver(double clf) : clf(clf) {}
    double minT;
    bool updateDomains(std::vector<Domain> &domains)
    {
        if (!fetchTimeStep(domains))
            return false;
        if (!updateX(domains))
            return false;
        if (!updateY(domains))
            return false;
        if (!updateZ(domains))
            return false;
        if (!updateZ(domains))
            return false;
        if (!updateY(domains))
            return false;
        if (!updateX(domains))
            return false;
        return true;
    }
    bool fetchTimeStep(std::vector<Domain> &domains)
    {
        minT = 1e10;
        for (auto &d : domains)
        {
            if (!timeStep(d))
                return false;
        }
        return true;
    }
    bool timeStep(Domain &d)
    {
        for (int i = 1; i < d.nx - 1; i++)
        {
            for (int j = 1; j < d.ny - 1; j++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    if (clf * d.boxDims / (d.u(i, j, k) + d.a(i, j, k)) < minT)
                        minT = clf * d.boxDims / (d.u(i, j, k) + d.a(i, j, k));
                    if (clf * d.boxDims / (d.v(i, j, k) + d.a(i, j, k)) < minT)
                        minT = clf * d.boxDims / (d.v(i, j, k) + d.a(i, j, k));
                    if (clf * d.boxDims / (d.w(i, j, k) + d.a(i, j, k)) < minT)
                        minT = clf * d.boxDims / (d.w(i, j, k) + d.a(i, j, k));
                }
            }
        }
        return true;
    }
    bool updateX(std::vector<Domain> &domains)
    {
        for (auto &d : domains)
        {
            if (!xFaces(d))
                return false;
            if (!xBoxes(d))
                return false;
            if (!xGhosts(d))
                return false;
        }
        return true;
    }
    bool updateY(std::vector<Domain> &domains)
    {
        for (auto &d : domains)
        {
            if (!yFaces(d))
                return false;
            if (!yBoxes(d))
                return false;
            if (!yGhosts(d))
                return false;
        }
        return true;
    }
    bool updateZ(std::vector<Domain> &domains)
    {
        for (auto &d : domains)
        {
            if (!zFaces(d))
                return false;
            if (!zBoxes(d))
                return false;
            if (!zGhosts(d))
                return false;
        }
        return true;
    }
    bool xFaces(Domain &d)
    {
        for (int i = 0; i < d.nxFaces; i++)
        {
            for (int j = 0; j < d.ny; j++)
            {
                for (int k = 0; k < d.nz; k++)
                {
                    if (!rs.findStar(d.rho(i, j, k), d.u(i, j, k), d.v(i, j, k), d.w(i, j, k), d.a(i, j, k), d.p(i, j, k), d.rho(i + 1, j, k), d.u(i + 1, j, k), d.v(i + 1, j, k), d.w(i + 1, j, k), d.a(i + 1, j, k), d.p(i + 1, j, k), d.xfAt(i, j, k)))
                        return false;
                }
            }
        }
        return true;
    }
    bool yFaces(Domain &d)
    {
        for (int i = 0; i < d.nx; i++)
        {
            for (int j = 0; j < d.nyFaces; j++)
            {
                for (int k = 0; k < d.nz; k++)
                {
                    if (!rs.findStar(d.rho(i, j, k), d.v(i, j, k), d.u(i, j, k), d.w(i, j, k), d.a(i, j, k), d.p(i, j, k), d.rho(i, j + 1, k), d.v(i, j + 1, k), d.u(i, j + 1, k), d.w(i, j + 1, k), d.a(i, j + 1, k), d.p(i, j + 1, k), d.yfAt(i, j, k)))
                        return false;
                }
            }
        }
        return true;
    }
    bool zFaces(Domain &d)
    {
        for (int i = 0; i < d.nx; i++)
        {
            for (int j = 0; j < d.ny; j++)
            {
                for (int k = 0; k < d.nzFaces; k++)
                {
                    if (!rs.findStar(d.rho(i, j, k), d.w(i, j, k), d.v(i, j, k), d.u(i, j, k), d.a(i, j, k), d.p(i, j, k), d.rho(i, j, k + 1), d.w(i, j, k + 1), d.v(i, j, k + 1), d.u(i, j, k + 1), d.a(i, j, k + 1), d.p(i, j, k + 1), d.zfAt(i, j, k)))
                        return false;
                }
            }
        }
        return true;
    }
    bool xBoxes(Domain &d)
    {
        for (int i = 1; i < d.nx - 1; i++)
        {
            for (int j = 1; j < d.ny - 1; j++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    BoxConserved b(d.rho(i, j, k), d.u(i, j, k), d.v(i, j, k), d.w(i, j, k), d.p(i, j, k), d.a(i, j, k));
                    double u1 = b.u1() + minT * (d.xfAt(i - 1, j, k).f1 - d.xfAt(i, j, k).f1) / d.boxDims;
                    double u2 = b.u2() + minT * (d.xfAt(i - 1, j, k).f2 - d.xfAt(i, j, k).f2) / d.boxDims;
                    double u3 = b.u3() + minT * (d.xfAt(i - 1, j, k).f3 - d.xfAt(i, j, k).f3) / d.boxDims;
                    double u4 = b.u4() + minT * (d.xfAt(i - 1, j, k).f4 - d.xfAt(i, j, k).f4) / d.boxDims;
                    double u5 = b.u5() + minT * (d.xfAt(i - 1, j, k).f5 - d.xfAt(i, j, k).f5) / d.boxDims;
                    b.updateFromConservatives(u1, u2, u3, u4, u5);
                }
            }
        }
        return true;
    }
    bool yBoxes(Domain &d)
    {
        for (int i = 1; i < d.nx - 1; i++)
        {
            for (int j = 1; j < d.ny - 1; j++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    BoxConserved b(d.rho(i, j, k), d.u(i, j, k), d.v(i, j, k), d.w(i, j, k), d.p(i, j, k), d.a(i, j, k));
                    double u1 = b.u1() + minT * (d.yfAt(i, j - 1, k).f1 - d.yfAt(i, j, k).f1) / d.boxDims;
                    double u2 = b.u2() + minT * (d.yfAt(i, j - 1, k).f2 - d.yfAt(i, j, k).f2) / d.boxDims;
                    double u3 = b.u3() + minT * (d.yfAt(i, j - 1, k).f3 - d.yfAt(i, j, k).f3) / d.boxDims;
                    double u4 = b.u4() + minT * (d.yfAt(i, j - 1, k).f4 - d.yfAt(i, j, k).f4) / d.boxDims;
                    double u5 = b.u5() + minT * (d.yfAt(i, j - 1, k).f5 - d.yfAt(i, j, k).f5) / d.boxDims;
                    b.updateFromConservatives(u1, u2, u3, u4, u5);
                }
            }
        }
        return true;
    }
    bool zBoxes(Domain &d)
    {
        for (int i = 1; i < d.nx - 1; i++)
        {
            for (int j = 1; j < d.ny - 1; j++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    BoxConserved b(d.rho(i, j, k), d.u(i, j, k), d.v(i, j, k), d.w(i, j, k), d.p(i, j, k), d.a(i, j, k));
                    double u1 = b.u1() + minT * (d.zfAt(i, j, k - 1).f1 - d.zfAt(i, j, k).f1) / d.boxDims;
                    double u2 = b.u2() + minT * (d.zfAt(i, j, k - 1).f2 - d.zfAt(i, j, k).f2) / d.boxDims;
                    double u3 = b.u3() + minT * (d.zfAt(i, j, k - 1).f3 - d.zfAt(i, j, k).f3) / d.boxDims;
                    double u4 = b.u4() + minT * (d.zfAt(i, j, k - 1).f4 - d.zfAt(i, j, k).f4) / d.boxDims;
                    double u5 = b.u5() + minT * (d.zfAt(i, j, k - 1).f5 - d.zfAt(i, j, k).f5) / d.boxDims;
                    b.updateFromConservatives(u1, u2, u3, u4, u5);
                }
            }
        }
        return true;
    }
    bool xGhosts(Domain &d)
    {
        // side 0 is the +x side
        if (d.sides[0]) // transmissive boundary condition
        {
            for (int j = 1; j < d.ny - 1; j++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    d.at(d.nx - 1, j, k) = d.sides[0]->at(0, j, k);
                }
            }
        }
        else // reflective boundary condition
        {
            for (int j = 1; j < d.ny - 1; j++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    d.at(d.nx - 1, j, k) = d.at(d.nx - 2, j, k);
                    d.at(d.nx - 1, j, k).u = -d.at(d.nx - 1, j, k).u;
                }
            }
        }
        // side 1 is the -x side
        if (d.sides[1]) // transmissive boundary condition
        {
            for (int j = 1; j < d.ny - 1; j++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    d.at(0, j, k) = d.sides[1]->at(d.sides[1]->nx - 1, j, k);
                }
            }
        }
        else // reflective boundary condition
        {
            for (int j = 1; j < d.ny - 1; j++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    d.at(0, j, k) = d.at(1, j, k);
                    d.at(0, j, k).u = -d.at(0, j, k).u;
                }
            }
        }
        return true;
    }
    bool yGhosts(Domain &d)
    {
        // side 2 is the +y side
        if (d.sides[2]) // transmissive boundary condition
        {
            for (int i = 1; i < d.nx - 1; i++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    d.at(i, d.ny - 1, k) = d.sides[2]->at(i, 0, k);
                }
            }
        }
        else // reflective boundary condition
        {
            for (int i = 1; i < d.nx - 1; i++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    d.at(i, d.ny - 1, k) = d.at(i, d.ny - 2, k);
                    d.at(i, d.ny - 1, k).v = -d.at(i, d.ny - 1, k).v;
                }
            }
        }
        // side 3 is the -y side
        if (d.sides[3]) // transmissive boundary condition
        {
            for (int i = 1; i < d.nx - 1; i++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    d.at(i, 0, k) = d.sides[3]->at(i, d.sides[3]->ny - 1, k);
                }
            }
        }
        else // reflective boundary condition
        {
            for (int i = 1; i < d.nx - 1; i++)
            {
                for (int k = 1; k < d.nz - 1; k++)
                {
                    d.at(i, 0, k) = d.at(i, 1, k);
                    d.at(i, 0, k).v = -d.at(i, 0, k).v;
                }
            }
        }
        return true;
    }
    bool zGhosts(Domain &d)
    {
        // side 4 is the +z side
        if (d.sides[4]) // transmissive boundary condition
        {
            for (int i = 1; i < d.nx - 1; i++)
            {
                for (int j = 1; j < d.ny - 1; j++)
                {
                    d.at(i, j, d.nz - 1) = d.sides[4]->at(i, j, 0);
                }
            }
        }
        else // reflective boundary condition
        {
            for (int i = 1; i < d.nx - 1; i++)
            {
                for (int j = 1; j < d.ny - 1; j++)
                {
                    d.at(i, j, d.nz - 1) = d.at(i, j, d.nz - 2);
                    d.at(i, j, d.nz - 1).w = -d.at(i, j, d.nz - 1).w;
                }
            }
        }
        // side 5 is the -z side
        if (d.sides[5]) // transmissive boundary condition
        {
            for (int i = 1; i < d.nx - 1; i++)
            {
                for (int j = 1; j < d.ny - 1; j++)
                {
                    d.at(i, j, 0) = d.sides[5]->at(i, j, d.sides[5]->nz - 1);
                }
            }
        }
        else // reflective boundary condition
        {
            for (int i = 1; i < d.nx - 1; i++)
            {
                for (int j = 1; j < d.ny - 1; j++)
                {
                    d.at(i, j, 0) = d.at(i, j, 1);
                    d.at(i, j, 0).w = -d.at(i, j, 0).w;
                }
            }
        }
        return true;
    }

private:
    RiemannSolver rs;
    double clf;
};

class FileHandler
{
public:
    FileHandler(std::string filename) : filename(filename)
    {
        file = H5::H5File(filename + ".h5", H5F_ACC_TRUNC);
    }
    std::string filename;
    H5::H5File file;
    void writeTimestep(std::vector<Domain> &domains, const double &time, const int &iteration)
    {
        std::string groupname = "timestep_" + std::to_string(iteration);
        H5::Group timestepGroup = file.createGroup(groupname.c_str());

        // write the time
        hsize_t dims[1] = {1};
        H5::DataSpace dataspace = H5::DataSpace(1, dims);
        H5::DataSet timeDataset = timestepGroup.createDataSet("time", H5::PredType::NATIVE_DOUBLE, dataspace);
        timeDataset.write(&time, H5::PredType::NATIVE_DOUBLE);

        // write the domain data
        for (int i = 0; i < domains.size(); i++)
        {
            std::string domainGroupname = "domain_" + std::to_string(i);
            H5::Group domainGroup = timestepGroup.createGroup(domainGroupname.c_str());
            writeDomain(domains[i], domainGroup);
        }
    }
    void writeDomain(Domain &d, H5::Group &domainGroup)
    {
        // write the domain dimensions
        hsize_t dims[3] = {d.nx, d.ny, d.nz};
        H5::DataSpace dataspace = H5::DataSpace(3, dims);
        H5::DataSet dataset = domainGroup.createDataSet("rho", H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(d.rho_.data(), H5::PredType::NATIVE_DOUBLE);
        dataspace = H5::DataSpace(3, dims);
        dataset = domainGroup.createDataSet("u", H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(d.u_.data(), H5::PredType::NATIVE_DOUBLE);
        dataspace = H5::DataSpace(3, dims);
        dataset = domainGroup.createDataSet("v", H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(d.v_.data(), H5::PredType::NATIVE_DOUBLE);
        dataspace = H5::DataSpace(3, dims);
        dataset = domainGroup.createDataSet("w", H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(d.w_.data(), H5::PredType::NATIVE_DOUBLE);
        dataspace = H5::DataSpace(3, dims);
        dataset = domainGroup.createDataSet("p", H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(d.p_.data(), H5::PredType::NATIVE_DOUBLE);
        dataspace = H5::DataSpace(3, dims);
        dataset = domainGroup.createDataSet("ghostCellMask", H5::PredType::NATIVE_UINT8, dataspace);
        dataset.write(d.ghostCellMask.data(), H5::PredType::NATIVE_UINT8);
    }
};

class XDMFHandler
{
public:
    XDMFHandler(std::vector<Domain> &domains, std::string filename, std::vector<double> &time) : filename(filename)
    {
        file.open(filename + ".xmf");
        writeHeader();
        writeGrids(domains, time);
        writeFooter();
        file.close();
    }

private:
    std::string filename;
    std::ofstream file;
    void writeHeader()
    {
        file << "<?xml version=\"1.0\" ?>\n";
        file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        file << "<Xdmf Version=\"2.0\">\n";
        file << "  <Domain>\n";
    }
    void writeGrids(std::vector<Domain> &domains, std::vector<double> &time)
    {
        file << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";
        for (int i = 0; i < time.size(); i++)
        {
            file << "      <Grid Name=\"Timestep_" << i << "\" GridType=\"Collection\" CollectionType=\"Spatial\">\n";
            file << "        <Time Value=\"" << time[i] << "\"/>\n";
            writeTimestep(domains, i);
            file << "      </Grid>\n";
        }
        file << "    </Grid>\n";
    }
    void writeFooter()
    {
        file << "  </Domain>\n";
        file << "</Xdmf>\n";
    }
    void writeTimestep(std::vector<Domain> &domains, const int &it)
    {
        int i = 0;
        for (auto &d : domains)
        {
            file << "        <Grid Name=\"Domain_" << i << "\" GridType=\"Uniform\">\n";
            file << "          <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"" << d.nx + 1 << " " << d.ny + 1 << " " << d.nz + 1 << "\"/>\n";
            file << "          <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
            file << "            <DataItem Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n";
            if (i == 0)
                file << "              0 0 0\n";
            else
                file << "              " << 2 << " 0 0\n"; // temp
            file << "            </DataItem>\n";
            file << "            <DataItem Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n";
            file << "              " << d.boxDims << " " << d.boxDims << " " << d.boxDims << "\n";
            file << "            </DataItem>\n";
            file << "          </Geometry>\n";
            file << "          <Attribute Name=\"rho\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
            file << "            <DataItem Format=\"HDF\" Dimensions=\"" << d.nx << " " << d.ny << " " << d.nz << "\" NumberType=\"Float\" Precision=\"8\">" << filename << ".h5:/timestep_" << it << "/domain_" << i << "/rho</DataItem>\n";
            file << "          </Attribute>\n";
            file << "          <Attribute Name=\"u\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
            file << "            <DataItem Format=\"HDF\" Dimensions=\"" << d.nx << " " << d.ny << " " << d.nz << "\" NumberType=\"Float\" Precision=\"8\">" << filename << ".h5:/timestep_" << it << "/domain_" << i << "/u</DataItem>\n";
            file << "          </Attribute>\n";
            file << "          <Attribute Name=\"v\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
            file << "            <DataItem Format=\"HDF\" Dimensions=\"" << d.nx << " " << d.ny << " " << d.nz << "\" NumberType=\"Float\" Precision=\"8\">" << filename << ".h5:/timestep_" << it << "/domain_" << i << "/v</DataItem>\n";
            file << "          </Attribute>\n";
            file << "          <Attribute Name=\"w\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
            file << "            <DataItem Format=\"HDF\" Dimensions=\"" << d.nx << " " << d.ny << " " << d.nz << "\" NumberType=\"Float\" Precision=\"8\">" << filename << ".h5:/timestep_" << it << "/domain_" << i << "/w</DataItem>\n";
            file << "          </Attribute>\n";
            file << "          <Attribute Name=\"p\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
            file << "            <DataItem Format=\"HDF\" Dimensions=\"" << d.nx << " " << d.ny << " " << d.nz << "\" NumberType=\"Float\" Precision=\"8\">" << filename << ".h5:/timestep_" << it << "/domain_" << i << "/p</DataItem>\n";
            file << "          </Attribute>\n";
            file << "          <Attribute Name=\"ghostCellMask\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
            file << "            <DataItem Format=\"HDF\" Dimensions=\"" << d.nx << " " << d.ny << " " << d.nz << "\" NumberType=\"UInt8\" Precision=\"8\">" << filename << ".h5:/timestep_" << it << "/domain_" << i << "/ghostCellMask</DataItem>\n";
            file << "          </Attribute>\n";
            file << "        </Grid>\n";
            i++;
        }
    }
};

int main()
{
    std::vector<Domain> domains(2);
    State s(1.0, 0.0, 0.0, 0.0, 1.0);
    domains[0].setup(1, 1, 1, 5, s);
    State s2(0.125, 0.0, 0.0, 0.0, 0.1);
    domains[1].setup(1, 1, 1, 5, s2);
    domains[0].sides[0] = &domains[1];
    domains[1].sides[1] = &domains[0];
    DomainSolver ds(0.5);
    std::vector<double> t(1, 0.0);
    double tEnd = 50;
    int iteration = 0;

    std::string filename("test");

    FileHandler fh(filename);
    fh.writeTimestep(domains, t.back(), iteration);

    while (t.back() < tEnd)
    {
        if (!(ds.updateDomains(domains)))
            break;
        t.push_back(t.back() + 2 * ds.minT);
        iteration++;
        std::cout << "Time elapsed: " << t.back() << std::endl;
        fh.writeTimestep(domains, t.back(), iteration);
    }
    XDMFHandler xh(domains, filename, t);
    std::cout << "done" << std::endl;
    std::cout << "pressure centre domain 0: " << domains[0].p(2, 2, 2) << std::endl;
    std::cout << "pressure centre domain 1: " << domains[1].p(2, 2, 2) << std::endl;
    return 0;
}