#ifndef SOLVERIEMANN_H
#define SOLVERIEMANN_H

#include <vector>

class SolveRiemann
{
public:
    static std::vector<double> getStars(double &pL, double &rhoL, double &uL, double &vL, double &aL, double &pR, double &rhoR, double &uR, double &vR, double &aR);

private:
    static double guessp(double &pL, double &rhoL, double &uL, double &vL, double &aL, double &pR, double &rhoR, double &uR, double &vR, double &aR);
    static std::vector<double> prefun(double &P, double &PK, double &CK, double &DK);
    static double getRho(double &u, double &p, double &pL, double &rhoL, double &pR, double &rhoR);
    static constexpr double gammma = 1.4;
    static constexpr double G1 = (gammma - 1) / (2 * gammma);
    static constexpr double G2 = (gammma + 1) / (2 * gammma);
    static constexpr double G3 = 2 * gammma / (gammma - 1);
    static constexpr double G4 = 2 / (gammma - 1);
    static constexpr double G5 = 2 / (gammma + 1);
    static constexpr double G6 = (gammma - 1) / (gammma + 1);
    static constexpr double G7 = (gammma - 1) / 2;
    static constexpr double G8 = gammma - 1;
};

#endif
