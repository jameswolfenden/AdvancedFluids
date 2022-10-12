#ifndef POINT_H
#define POINT_H

class Point
{
public:
    double p, rho, u, a;
    Point() {}
    double aCalc();
    Point(double p, double rho, double u);
};


#endif
