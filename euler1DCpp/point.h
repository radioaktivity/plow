#ifndef point_H
#define point_H

#include <iostream>
#include <vector>
#include <cmath>

class Point{
    public:
        float X;
        float Y;
        Point();
        Point(float x, float y);
        void set_coordinates(float x, float y);
        std::vector<float> getVecBetween(Point p_other);
        float distance(Point p_other);
};

#endif // point_H
