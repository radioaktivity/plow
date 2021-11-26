#include "point.h"
#include "face.h"
#include "cell.h"
#include "convert.h"

#include <iostream>
#include <vector>


using namespace std;

Cell::Cell(){
  // Do nothing
}

Cell::Cell(int NUMBER, std::vector<Point*> BOUNDARY_POINTS){
  number = NUMBER;

  // setting the boundary points 
  for (Point* p : BOUNDARY_POINTS)
  {
    boundary_points.push_back(p);
  }
  Cell::calc_all();

}

void Cell::calc_all(){
  Cell::calc_center();
  Cell::calc_length();
}

void Cell::calc_length(){
  length = boundary_points.operator[](0)->distance(*boundary_points.operator[](1));
  volume = length;
}

void Cell::calc_center(){
  float x = 1./2. * (boundary_points.operator[](0)->X +  boundary_points.operator[](1)->X);
  float y = 1./2. * (boundary_points.operator[](0)->Y +  boundary_points.operator[](1)->Y);
  center.set_coordinates(x,y);
}

void Cell::calc_primitive(){
  float * result;
  result = getPrimitive(Mass, Momx, Energy, volume);
  rho = result[0];
  u = result[1];
  p = result[2];
}
void Cell::calc_conserved(){
  float * result;
  result =  getConserved(rho, u, p, volume);
  Mass = result[0];
  Momx = result[1];
  Energy = result[2];
}

void Cell::calc_gradients(string type){

  if (type == "central"){
    float d = 2* length;
    rho_dx = (neighbors.operator[](1)->rho - neighbors.operator[](0)->rho)/d;
    u_dx = (neighbors.operator[](1)->u - neighbors.operator[](0)->u)/d;
    p_dx = (neighbors.operator[](1)->p - neighbors.operator[](0)->p)/d;
  }
  else{
    throw "No other gradients implemented";
  }

}

void Cell::extrapol_in_time(float dt){
  float gamma = 5./3.;
  rho = rho - 0.5* dt * (u * rho_dx + rho * u_dx);
  u = u - 0.5 * dt * (u * u_dx + 1/rho * p_dx);
  p = p - 0.5 * dt * (gamma * p * u_dx + u * p_dx);
}

void Cell::extrapol_2_faces(){
  // left Face
  Face* face_L = faces.operator[](0);
  // right Face
  Face* face_R = faces.operator[](1);
  std::vector<float> fn = center.getVecBetween(face_L->center);
  float rho_face = rho + rho_dx * fn.operator[](0);
  float u_face = u + u_dx * fn.operator[](0);
  float p_face = p + p_dx * fn.operator[](0);

  // give the left face is values
  // on pov of left face cell is on side = "R"
  face_L->getPrimitives(rho_face,u_face,p_face, 1);

  rho_face = rho - rho_dx * fn.operator[](0);
  u_face = u - u_dx * fn.operator[](0);
  p_face = p - p_dx * fn.operator[](0);
  face_R->getPrimitives(rho_face,u_face,p_face, 0);
}

void Cell::apply_flux(float dt){
  Face* face_L = faces.operator[](0);
  Face* face_R = faces.operator[](1);
  Mass -= dt * (face_L->flux_Mass-face_R->flux_Mass) * volume;
  Momx -= dt * (face_L->flux_Momx-face_R->flux_Momx) * volume;
  Energy -= dt * (face_L->flux_Energy-face_R->flux_Energy) * volume;
  Cell::calc_primitive();
}

void Cell::print_cell(){
  std::cout<< " Cell number " << number << " Center x: " << center.X << std::endl;
}

void Cell::print_neighbors(){
    std::cout<< " Cell number " << number << " neighbors L, R: " << neighbors.operator[](0)->number << ", " << neighbors.operator[](1)->number << std::endl;
}
void Cell::print_faces(){
  std::cout<< " Cell number " << number << " Faces L, R: " << faces.operator[](0)->number << ", "  << faces.operator[](1)->number << " Face L Center "<< faces.operator[](0)->center.X << std::endl;
}


int cell_test() {
    Point p1(0.0, 0.0);
    Point p2(0.5, 0.0);
    Point p3(1.0, 0.0);
    std::vector<Point*> vec1 = {&p1, &p2};
    std::vector<Point*> vec2 = {&p2, &p3};

    Cell c1(0, vec1);
    Cell c2(0, vec2);
    c1.calc_center();
    c2.calc_center();
    std::cout << c1.boundary_points.operator[](0)->X << std::endl;

    std::cout << " Center of Cell 1 x: " << c1.center.X << " y: " << c1.center.Y << std::endl;
    std::cout << " Center of Cell 2 x: " << c2.center.X << " y: " << c2.center.Y << std::endl;


    return 0;
}
