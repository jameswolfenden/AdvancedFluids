#include <cmath>
#include <iostream>
#include "Domain.h"
#include "../fluidConsts.h"

Domain::Domain(double xPhysical, int pointsSize)
{
    this->xPhysical = xPhysical; // set physical size
    this->pointsSize = pointsSize;
    halfSize = pointsSize - 1;
    xBox = xPhysical / pointsSize; // work out size of each box
    elapsedTime = 0;
    points.resize(pointsSize);
    halfPoints.resize(halfSize);

    for (int i = 1; i < (pointsSize - 1); i++) // exclude ghost cells
        points[i].updatePrimatives(1, 1, 0);

    for (int i = 1; i < (pointsSize / 2) + 1; i++) // set first half of cells
        points[i].updatePrimatives(0.1, 0.125, 0);

    points[0].updatePrimatives(points[1].p, points[1].rho, -points[1].u);                                                     // ghost points
    points[pointsSize - 1].updatePrimatives(points[pointsSize - 2].p, points[pointsSize - 2].rho, -points[pointsSize - 2].u); // ghost points
}
void Domain::findHalfs()
{
    for (int i = 0; i < halfSize; i++)
    {
        // each side is a pointer to the cell on each side of the half point
        Point *sides[2];
        sides[0] = &points[i];
        sides[1] = &points[i + 1];

        halfPoints[i].findStar(sides); // find the star values for the half point

        if (domainOutput)
            std::cout << "start: " << points[i].p << " mid: " << halfPoints[i].p << " end: " << points[i + 1].p << std::endl;
    }
}
void Domain::updatePoints()
{
    minT = timeStep();
    findHalfs();
    for (int i = 1; i < (pointsSize - 1); i++)
    {
        // find conservatives from the half point fluxes
        double u1 = points[i].u1() + minT / xBox * (halfPoints[i - 1].f1() - halfPoints[i].f1());
        double u2 = points[i].u2() + minT / xBox * (halfPoints[i - 1].f2() - halfPoints[i].f2());
        double u3 = points[i].u3() + minT / xBox * (halfPoints[i - 1].f3() - halfPoints[i].f3());
        points[i].updateConservatives(u1, u2, u3); // update the primatives in the points array from the conservatives found
        if (domainOutput)
            std::cout << "new p : " << points[i].p << " made of : " << halfPoints[i - 1].p << " and: " << halfPoints[i].p << std::endl;
    }
    points[0].updatePrimatives(points[1].p, points[1].rho, -points[1].u);                                                     // ghost points
    points[pointsSize - 1].updatePrimatives(points[pointsSize - 2].p, points[pointsSize - 2].rho, -points[pointsSize - 2].u); // ghost points
    elapsedTime += minT;
}

double Domain::timeStep()
{
    double step = 1000000; // high value so its larger then all starting, should change to some inital calc
    for (int i = 1; i < (pointsSize - 1); i++)
    {
        double C_clf = 0.49;
        // this needs optimisation
        if (C_clf * xBox / (std::abs(points[i].u) + points[i].a) < step) // set the deltaT for the next iteration being 0.9 * the min time step for this iteration
            step = C_clf * xBox / (std::abs(points[i].u) + points[i].a);
    }
    return step;
}