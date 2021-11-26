#ifndef face_H
#define face_H

#include "point.h"

#include <iostream>
#include <vector>
#include <cmath>

class Face{
    public:
    int number;
    Face* wormhole_face;
    bool isL;
    bool isR;
    
    Point center;
    float rho_L;
    float rho_R;
    float u_L;
    float u_R;
    float p_L;
    float p_R;

    float flux_Mass;
    float flux_Momx;
    float flux_Energy;
    
    Face();
    Face(int NUMBER, Point center);
    void getPrimitives(float rho, float u, float p, int side);
    void calcFlux();
    void getFlux();
};

#endif
