#include "create_mesh.h"
#include <iostream>

using namespace std;

void create_mesh(int n, Cell* cells, Face * faces){

  float dimension = 1.0;

  for ( int i=0; i<n; i++){
    Point p1(((float)i) /n * dimension, 0.0);
    Point p2(((float)i+1) /n * dimension, 0.0);
    std::vector<Point*> vec = {&p1, &p2};
    Cell c_new(i, vec); 
    cells[i] = c_new;     
    Face f_new(i, p1);
    faces[i] = f_new;
  }
  faces[n] = Face(n,Point(dimension, 0.0));
  faces[0].wormhole_face = &faces[n];
  faces[n].wormhole_face = &faces[0];

  for ( int i=0; i<n; i++){
    Cell c = cells[i];
    c.faces.push_back(&faces[i]);
    c.faces.push_back(&faces[i+1]);
  }

  for ( int i=1 ; i<n-1; i++){
    // Add left neighbor
    cells[i].neighbors.push_back(&cells[i-1]);
    // Add right neighbor
    cells[i].neighbors.push_back(&cells[i+1]);
  }

  cells[0].neighbors.push_back(&cells[n-1]);
  cells[0].neighbors.push_back(&cells[1]);
  cells[n-1].neighbors.push_back(&cells[n-2]);
  cells[n-1].neighbors.push_back(&cells[0]);

  /*
  for ( int i=0 ; i<n; i++){
    cells[i].calc_all();
  }
  */

}