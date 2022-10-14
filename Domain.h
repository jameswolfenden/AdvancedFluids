#ifndef DOMAIN_H
#define DOMAIN_H

#include "WavesDataPoint.h"
#include <vector>

const int pointsSize = 1000;
const int halfSize = pointsSize-1;

class Domain{
    public:
    int xPhysical;
    double xBox;
    double minT;
    double elapsedTime;
    std::vector<Point> points;
    std::vector<WavesDataPoint> halfPoints;
    Domain(double xPhysical);
    void updatePoints();
    void findHalfs();
};


#endif

