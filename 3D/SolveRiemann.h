#ifndef SOLVERIEMANN_H
#define SOLVERIEMANN_H

class SolveRiemann{
    public:
    bool findStar(const double &rhoL, const double &uL, const double &vL, const double &wL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &wR, const double &aR, const double &pR, double &rho, double &u, double &v, double &w, double &a, double &p);
    private:
    bool testVacuum(const double &rhoL, const double &uL, const double &vL, const double &wL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &wR, const double &aR, const double &pR, double &rho, double &u, double &v, double &w, double &a, double &p);
    bool pickStartVal(const int errorStage, const double &rhoL, const double &uL, const double &vL, const double &wL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &wR, const double &aR, const double &pR, double &rho, double &u, double &v, double &w, double &a, double &p);
    void iterateP(const double &rhoL, const double &uL, const double &vL, const double &wL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &wR, const double &aR, const double &pR, double &rho, double &u, double &v, double &w, double &a, double &p);
};


#endif

