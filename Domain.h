#ifndef DOMAIN_H
#define DOMAIN_H

#include "WavesDataPoint.h"

const int pointsSize = 32;
const int halfSize = pointsSize-1;

class Domain{
    public:
    int xPhysical;
    double xBox;
    double minT;
    Point points[pointsSize];
    WavesDataPoint halfPoints[halfSize];
    Domain(double xPhysical);
    void updatePoints();
    void findHalfs();
};


#endif

