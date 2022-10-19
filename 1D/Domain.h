#ifndef DOMAIN_H
#define DOMAIN_H

#include "Point.h"
#include <vector>

class Domain{
    public:
    int xPhysical;
    double xBox;
    double minT;
    double elapsedTime;
    int pointsSize, halfSize;
    std::vector<Point> points;
    std::vector<Point> halfPoints;
    Domain(double xPhysical, int pointsSize);
    void updatePoints();
    void findHalfs();
};


#endif

