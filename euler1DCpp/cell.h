#ifndef cell_H
#define cell_H

#include <iostream>
#include <vector>
#include <cmath>


using namespace std;
// Create a Car class with some attributes
class Cell {
  public:
    // Organisation
    int number;   
    std::vector<Point*> boundary_points;
    std::vector<Cell*> neighbors;
    std::vector<Face*> faces;

    // Geometrie
    Point center;
    float volume;
    float length;

    // primitives
    float rho;
    float u; 
    float p;

    // gradients
    float rho_dx;
    float u_dx;
    float p_dx;

    // conservatives
    float Mass;
    float Momx;
    float Energy;
    
    Cell();
    Cell(int NUMBER, std::vector<Point*> BOUNDARY_POINTS);
    
    void calc_primitive();
    void calc_conserved();
    
    void calc_center();
    void calc_length();
    void calc_all();

    void calc_gradients(string type);

    void extrapol_in_time(float dt);

    void extrapol_2_faces();

    void apply_flux(float dt);

    void print_cell();

    void print_neighbors();

    void print_faces();


};

int cell_test();

#endif