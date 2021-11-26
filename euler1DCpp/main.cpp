
#include "algos.h"
#include "convert.h"

#include "point.h"
#include "face.h"
#include "cell.h"
#include "create_mesh.h"


#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>

using namespace std;

int main(){

  float t = 0.;
  float t_end = 10.;
  float dt;
  int n = 20;
  float courant_fac = 0.2;
  float gamma = 5./3.;

  Cell cells[n];   
  Face faces[n+1];
  create_mesh(n, cells, faces);

  for(int i = 0; i<n; i++){
    Cell * c = &cells[i];
    c->print_cell();
    c->print_neighbors();
    c->print_faces();

    if(c->number < n-10){
      c->rho = 1.;
      c->u = 1.;
      c->p = 2.5;
    }
    else {
      c->rho = 1.;
      c->u = 0.;
      c->p = 2.5;   
    }
  }
   for(Cell c : cells){
    std::cout<<c.u<<std::endl;
    }

  Cell* c1 = &cells[0];

  float possible_dts[n];

  for (int i = 0; i<n; i++){
    Cell * c = &cells[i];
    possible_dts[i] = courant_fac * c->length / (sqrt(gamma * c->p/c->rho) + c->u * c->u) ;
  }

  dt = get_min_from_array(possible_dts, n);

  ofstream myfile;
  myfile.open("U.csv");

  // Write to file
  for(Cell c : cells){
    myfile << c.u << "; ";
  }
  myfile << "; " << dt;
  myfile << std::endl;
  

  while (t<t_end){


    for(int i = 0; i<n; i++){
      Cell * c = &cells[i];
      c->calc_gradients("central");
    }
    for(int i = 0; i<n; i++){
      Cell * c = &cells[i];
      c->extrapol_in_time(dt);
    }
    for(int i = 0; i<n; i++){
      Cell * c = &cells[i];
      c->extrapol_2_faces();
    }
    
    for(int i = 0; i<=n; i++){
      Face * f = &faces[i];
      f->calcFlux();
    }
    for(int i = 0; i<n; i++){
      Cell * c = &cells[i];
      c->apply_flux(dt);
    }

    // Write to file
    for(Cell c : cells){
      myfile << c.u << "; ";
    }
    myfile << "; " << dt;
    myfile << std::endl;
    

    // Calculate new Timestep
    for (int i = 0; i<n; i++){
      Cell * c = &cells[i];
      possible_dts[i] = courant_fac * c->length / (sqrt(gamma * c->p/c->rho) + c->u * c->u) ;
    }

    dt = get_min_from_array(possible_dts, n);  
    t += dt;
  }
  myfile.close();
  return 0;

}