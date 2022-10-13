#ifndef POINT_H
#define POINT_H

class Point
{
public:
    double p, rho, u, a;
    Point(){}
    double aCalc(), f1(), f2(), f3(), u1(), u2(), u3();
    Point(double p, double rho, double u);
    void updatePrimatives(double p, double rho, double u);
    void updateConservatives(double u1, double u2, double u3);
};

#endif
