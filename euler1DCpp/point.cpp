
#include "point.h"

using namespace std;

class Point;

Point::Point(){

}

Point::Point(float x, float y){
            X = x;
            Y = y;
        }

vector<float> Point::getVecBetween(Point p_other)
{
    vector<float> vec (2);
    vec.operator[](0) = p_other.X - X;
    vec.operator[](1) = p_other.Y - Y;
    return vec;
}

void Point::set_coordinates(float x, float y){
    X = x;
    Y = y;
}

float Point::distance(Point p_other)
{
    float dx =  p_other.X - X;
    float dy =  p_other.Y - Y;
    return sqrt(dx*dx+dy*dy);
}

/*
int test_point()
{
    Point p1(1., 0.);
    std::cout << p1.X << " " << p1.Y << std::endl;
    Point p2(0., 1.);
    std::cout << p2.X << " " << p2.Y << std::endl;

    std::cout << p1.getVecBetween(p2).operator[](0) << " "<< p1.getVecBetween(p2).operator[](1) << std::endl;
    std::cout << p1.distance(p2) << std::endl;
    

    return 0;
}
*/