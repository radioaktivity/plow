#ifndef create_mesh_H
#define create_mesh_H

#include "point.h"
#include "face.h"
#include "cell.h"
#include "convert.h"

#include <iostream>
#include <vector>



void create_mesh(int n, Cell * cells, Face * faces);

#endif