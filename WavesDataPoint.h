#ifndef WAVESDATAPOINT_H
#define WAVESDATAPOINT_H

#include "Point.h"

class WavesDataPoint : public Point
{
public:
    double uShock, uHead, uTail;
    double rhos[2];
    double rhoFan[100], uFan[100], pFan[100], us[100];

    void findStar(Point *sides[]);
    void waveData(Point *sides[]);
    using Point::Point;
};


#endif
