import csv
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

from global_proporties import *

#plt.rcParams.update({'font.size': 50})
# from scipy.ndimage.interpolation import zoom

from cell_new import *
from point import *
from face_new import *


def create_mesh_rect(n=8, plot_cells=False):

    # Dimensions of flowfield [m]
    dimension_x = atm.dimension_x
    dimension_y = atm.dimension_y

    # Create point map and point array
    x = np.linspace(0,dimension_x,n)
    y = np.linspace(0,dimension_y,n)
    X, Y = np.meshgrid(x,y)
    points = np.vstack([X.ravel(), Y.ravel()]).T

    points_obj = []

    for p in points:
        points_obj.append(Point(p))


    faces_set = set()
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

        c.create_faces()

    if True:
        n = n-1
        for i, c in enumerate(cells):
            # NOTE: Neighbors are sorted always in order LEFT RIGHT TOP BOTTOM

            # Add left neighbor
            if (i%(n))==0: # left border
                c.add_neighbor(cells[i+n-1])
                c.faces[0].wormhole_face = cells[i+n-1].faces[1]
            else:
                c.add_neighbor(cells[i-1])
                
            # Add right neighbor
            if ((i+1)%n == 0):
                c.add_neighbor(cells[i-(n-1)])
                c.faces[1].wormhole_face = cells[i-(n-1)].faces[0]
            else:
                c.add_neighbor(cells[i+1])

            # Add top neighbor
            try:
                c.add_neighbor(cells[i+n])
            except:
                c.add_neighbor(cells[i-((n-1)*n)])
                c.faces[2].wormhole_face = cells[i-((n-1)*n)].faces[3]
            

            c.add_neighbor(cells[i-n])
            if (i-n)<0:
                c.faces[3].wormhole_face = cells[i-n].faces[2]

            if False:
                c.add_neighbor(cells[i+((n-1)*n)])
                c.faces[3].wormhole_face = cells[i+((n-1)*n)].faces[2]

    i = 0
    for c in cells:
        c.face_neighbor_check()
        for f in c.faces:
            faces_set.add(f)
        c.calc_all()
        #print('_'*50)
        #print("# Cell: ",c.number)
        #print('neighbors: ', c.neighbors[0].number, c.neighbors[1].number, c.neighbors[2].number, c.neighbors[3].number)
        #print('Voluume: ',c.volume)

    for f in faces_set:
        f.number = i
        i += 1
    
    if plot_cells:
        
        plt.plot(points[:,0], points[:,1], 'o')
        for j, p in enumerate(points):

            plt.text(p[0], p[1], j, ha='right') # label the points

        for c in cells:
            plt.text(c.center.X, c.center.Y, '#%d' % c.number, ha='center') # label triangles

        for i,f in enumerate(faces_set):

            plt.text(f.center.X, f.center.Y, '#%d' % f.number, ha='center') # label triangles
            #print(c.neighbors[0].number, c.neighbors[1].number, c.neighbors[2].number, c.neighbors[3].number)
            f.plot_border()

    return [cells, points_obj, faces_set]

if __name__ == "__main__":
    cells, points_obj, faces_set = create_mesh_rect(n=9,plot_cells=True)
    
    for c in cells: 
        print('_'*50)
        print('Cell number: ', c.number)
        print('center', c.center)
        print('Neighbors ', [n.number for n in c.neighbors])
        print('Faces ',[[f.number, f.center.__str__()] for f in c.faces])


    for f in faces_set:
        print('_'*50)
        print('Face number: ', f.number)
        try:
            print('wormhole Face: ', f.wormhole_face.number)
        except:
            print('wormhole Face: ', f.wormhole_face)
        print('cells connected: ', [c.number for c in f.cells_connected])

    plt.show()