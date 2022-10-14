#include <cmath>
#include <iostream>
#include "Domain.h"
#include "fluidConsts.h"

Domain::Domain(double xPhysical)
{
    this->xPhysical = xPhysical;
    xBox=xPhysical/pointsSize;
    minT = 9999; //initial
    elapsedTime = 0;
    points.resize(pointsSize);
    halfPoints.resize(halfSize);

    for (int i = 1; i < (pointsSize-1); i++) // exclude ghost cells
    {
        points[i].updatePrimatives(1,1,0);
    }
    for (int i=1; i<(pointsSize/2)+1;i++)
        points[i].updatePrimatives(0.1,0.125,0);
    points[0].updatePrimatives(points[1].p, points[1].rho, -points[1].u); // ghost points
    points[pointsSize-1].updatePrimatives(points[pointsSize-2].p, points[pointsSize-2].rho, -points[pointsSize-2].u); // ghost points
}
void Domain::findHalfs()
{
    for (int i = 0; i < halfSize; i++)
    {
        Point *sides[2];
        sides[0] = &points[i];
        sides[1] = &points[i+1];
        halfPoints[i].findStar(sides);
        if (domainOutput)
            std::cout << "start: "<< points[i].p << " mid: " << halfPoints[i].p << " end: " << points[i+1].p << std::endl;
        if (0.5*xBox/(halfPoints[i].uShock+halfPoints[i].a)<minT)
            minT = 0.5*xBox/(halfPoints[i].uShock+halfPoints[i].a);
    }

}
void Domain::updatePoints()
{
    findHalfs();
    for (int i = 1; i < (pointsSize-1); i++)
    {
        double u1 = points[i].u1() + minT/xBox * (halfPoints[i-1].f1() - halfPoints[i].f1());
        double u2 = points[i].u2() + minT/xBox * (halfPoints[i-1].f2() - halfPoints[i].f2());
        double u3 = points[i].u3() + minT/xBox * (halfPoints[i-1].f3() - halfPoints[i].f3());
        points[i].updateConservatives(u1,u2,u3);
        elapsedTime += minT;
        if (domainOutput)
            std::cout << "new p : "<< points[i].p << " made of : " << halfPoints[i-1].p << " and: " << halfPoints[i].p << std::endl;
    }
    points[0].updatePrimatives(points[1].p, points[1].rho, -points[1].u); // ghost points
    points[pointsSize-1].updatePrimatives(points[pointsSize-2].p, points[pointsSize-2].rho, -points[pointsSize-2].u); // ghost points

}