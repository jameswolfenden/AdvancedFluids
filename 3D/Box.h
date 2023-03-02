#ifndef BOX_H
#define BOX_H

class Box
{
public:
    double p, rho, u, v, w, a;
    Box(){}
    double aCalc(), f1(), f2(), f3(), f4(), f5(), u1(), u2(), u3(), u4(), u5(), g1(), g2(), g3(), g4(), g5(), h1(), h2(), h3(), h4(), h5();
    Box(double p, double rho, double u, double v, double w);
    void updatePrimatives(double p, double rho, double u, double v, double w);
    void updateConservatives(double u1, double u2, double u3, double u4, double u5);
};

#endif
