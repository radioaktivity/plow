import csv
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import zoom

from cell import *
from point import *
from face import *


def create_mesh(nx=4,ny=4,write_mesh=False, plot_cells=False):

    #Your statements here

    # Create point map as well as delaunay triangles

    x = np.linspace(0,1,nx)
    y = np.linspace(0,1,ny)


    i = 0
    j = 0
    k = 0

    # Create field of points
    print(f"Creating Delaunay")
    print("-"*50)
    points = np.zeros([len(x)*len(y), 2])
    k = 0

    for i in x:
        for j in y:
            if (k % 2) == 0:
                corr = 0.1
            else:
                corr = 0
            points[k,0] = i + corr
            points[k,1] = j 
            k += 1
    # run delaunay algorithm
    tri = Delaunay(points)

    # Write cells and points in csv for later consumption
    if write_mesh:
        with open('MESH/eggs.csv', 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(['CellNumber', 
                'point1x','point1y', 
                'point2x','point2y',
                'point3x','point3y',
                'neighbor1', 'neighbor2','neighbor3'])
            i = 0
            for i in range(len(tri.simplices)):
                cis = tri.simplices[i] # CornerIndicieS
                neighbors = tri.neighbors[i]
                spamwriter.writerow([i, 
                            points[cis[0]][0],points[cis[0]][1],
                            points[cis[1]][0],points[cis[1]][1],
                            points[cis[2]][0],points[cis[2]][1],
                            neighbors[0], neighbors[1], neighbors[2]])
                i += 1

    # Convert triangles in cell-objects and points in point-objects
    # ----------------------------------------------------------

    print("Creating Cells and Points")
    print("-"*50)
    # Cretae as many empty cells as delauny triangles
    cells = []
    for i in range(len(tri.simplices)):
        cells.append(Cell(i))

    # create a Point() instance for every point there is 
    points_obj = []
    for i in range(len(points)):
        new_point = Point(points[i])
        points_obj.append(new_point)


    print("Setting boundary points and neighbors")
    print("-"*50)

    # give every cell his neighbors and boundary points 
    for i in range(len(tri.simplices)):

        c = cells[i] # current cell

        bps = [] # List of boundary points for current cell
        for index in tri.simplices[i]: # parse the point array of the delauny tri 
            bps.append(points_obj[index])
        c.set_boundary_points(bps)
        c.calc_center()

        for new_neighbor in tri.neighbors[i]: # parse the neighbor array of the delauny tri 
            if not(new_neighbor == -1): # when entry is -1, current face has no neighbor
                c.add_neighbor(cells[new_neighbor])


    print("Creating Faces and Neighbor check. Calculating vectors")
   
    print("-"*50)
    faces_set = set()
    for c in cells:
        c.create_faces()
        c.face_neighbor_check()
        c.calc_all()
        for f in c.faces:
            faces_set.add(f)


    if plot_cells:
        print("Plotting Cells")
        print("-"*50)
        plt.triplot(points[:,0], points[:,1], tri.simplices)

        plt.plot(points[:,0], points[:,1], 'o')
        for j, p in enumerate(points):

            plt.text(p[0], p[1], j, ha='right') # label the points

        for j, s in enumerate(tri.simplices):

            p = points[s].mean(axis=0)

            plt.text(p[0], p[1], '#%d' % j, ha='center') # label triangles

        plt.xlim(-0.5, 1.5); plt.ylim(-0.5, 1.5)


    return [cells, points_obj, faces_set]