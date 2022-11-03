#ifndef CELL_H
#define CELL_H

#include <vector>

class Cell
{
public:
    double p, rho, u, v, a;
    Cell(){}
    double aCalc(), f1(), f2(), f3(), f4(), u1(), u2(), u3(), u4(), g1(), g2(), g3(), g4();
    Cell(double p, double rho, double u, double v);
    void updatePrimatives(double p, double rho, double u, double v);
    void updatePrimatives(std::vector<double> results);
    void updateConservatives(double u1, double u2, double u3, double u4);
    void xFindStar(Cell *sides[]);
    void yFindStar(Cell *sides[]);
};

#endif
