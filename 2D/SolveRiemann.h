#ifndef SOLVERIEMANN_H
#define SOLVERIEMANN_H

class SolveRiemann{
    public:
    bool findStar(const double &rhoL, const double &uL, const double &vL, const double &aL, const double &pL, const double &dyeL, const double &rhoR, const double &uR, const double &vR, const double &aR, const double &pR, const double &dyeR, double &rho, double &u, double &v, double &a, double &p , double &dye);
    private:
    bool testVacuum(const double &rhoL, const double &uL, const double &vL, const double &aL, const double &pL, const double &dyeL, const double &rhoR, const double &uR, const double &vR, const double &aR, const double &pR, const double &dyeR, double &rho, double &u, double &v, double &a, double &p, double &dye);
    bool pickStartVal(const int errorStage, const double &rhoL, const double &uL, const double &vL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &aR, const double &pR, double &rho, double &u, double &v, double &a, double &p);
    void iterateP(const double &rhoL, const double &uL, const double &vL, const double &aL, const double &pL, const double &rhoR, const double &uR, const double &vR, const double &aR, const double &pR, double &rho, double &u, double &v, double &a, double &p);
};


#endif

