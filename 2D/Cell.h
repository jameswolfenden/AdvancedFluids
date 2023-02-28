#ifndef CELL_H
#define CELL_H

class Cell
{
public:
    double p, rho, u, v, a, dye;
    Cell(){}
    double aCalc(), f1(), f2(), f3(), f4(), f5(), u1(), u2(), u3(), u4(), u5(), g1(), g2(), g3(), g4(), g5();
    Cell(double p, double rho, double u, double v, double dye);
    void updatePrimatives(double p, double rho, double u, double v, double dye);
    void updateConservatives(double u1, double u2, double u3, double u4, double u5);
};

#endif
