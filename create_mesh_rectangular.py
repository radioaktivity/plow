import csv
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

from global_proporties import *

plt.rcParams.update({'font.size': 10})
from scipy.ndimage.interpolation import zoom

from cell_new import *
from point import *
from face_new import *


def create_mesh_rect(n=8, plot_cells=False):

    # Dimensions of flowfield [m]
    dimension_x = atm.dimension_x
    dimension_y = atm.dimension_y

    # Create point map and point array
    x = np.linspace(0,100,n)
    y = np.linspace(0,100,n)
    X, Y = np.meshgrid(x,y)
    points = np.vstack([X.ravel(), Y.ravel()]).T

    points_obj = []

    for p in points:
        points_obj.append(Point(p))


    number_of_cells = (n-1)**2
    cells = []
    j=-1
    for i in range(0,number_of_cells):
        if i%(n-1)==0:
            j+=1
        c = Cell(i)
        bps = []
        bps.append(points_obj[i+j])
        bps.append(points_obj[i+1+j])
        bps.append(points_obj[i+n+j])
        bps.append(points_obj[i+n+1+j])

        c.set_boundary_points(bps)
        c.calc_center()

        cells.append(c)

    if True:
        n = n-1
        for i, c in enumerate(cells):
            
            
            # NOTE: Neighbors are sorted always in order LET RIGHT TOP BOTTOM

            # Add left neighbor
            if (i%(n))==0: # left border
                c.add_neighbor(cells[i+n-1])
            else:
                c.add_neighbor(cells[i-1])

            # Add right neighbor
            if ((i+1)%n == 0):
                c.add_neighbor(cells[i-(n-1)])
            else:
                c.add_neighbor(cells[i+1])


            # Add top neighbor
            try:
                c.add_neighbor(cells[i+n])
            except:
                c.add_neighbor(cells[i-((n-1)*n)])
            
            # Add bottom neighbor
            try:
                c.add_neighbor(cells[i-n])
            except:
                c.add_neighbor(cells[i+((n-1)*n)])


    faces_set = set()
    for c in cells:
        c.create_faces()
    for c in cells:
        c.face_neighbor_check()
        for f in c.faces:
            faces_set.add(f)
        c.calc_all()
        #print('_'*50)
        #print("# Cell: ",c.number)
        #print('neighbors: ', c.neighbors[0].number, c.neighbors[1].number, c.neighbors[2].number, c.neighbors[3].number)
        #print('Voluume: ',c.volume)
    
    if plot_cells:
        plt.plot(points[:,0], points[:,1], 'o')
        for j, p in enumerate(points):

            plt.text(p[0], p[1], j, ha='right') # label the points

        for c in cells:
            plt.text(c.center.X, c.center.Y, '#%d' % c.number, ha='center') # label triangles

        for i,f in enumerate(faces_set):

            plt.text(f.center.X, f.center.Y, '#%d' % i, ha='center') # label triangles
            #print(c.neighbors[0].number, c.neighbors[1].number, c.neighbors[2].number, c.neighbors[3].number)
            f.plot_border()
        

    return [cells, points_obj, faces_set]

if __name__ == "__main__":
    create_mesh_rect()