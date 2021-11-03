import csv
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from cell import *
from point import *
from scipy.ndimage.interpolation import zoom

def create_mesh(nx=4,ny=4,write_mesh=False, plot_cells=False):
    # Create point map as well as delaunay triangles

    x = np.linspace(0,1,nx)
    y = np.linspace(0,1,ny)

    # For demonstration give every cell a random value of U 
    # arr = np.random.uniform(size=(nx,nx))
    # arr_u = zoom(arr, nx)
    # arr = np.random.uniform(size=(ny,ny))
    # arr_v = zoom(arr, ny)
    # velocity_vector = np.zeros([len(x)*len(y), 2])

    i = 0
    j = 0
    k = 0
    # for i in range(arr_u.shape[0]):
    #     for j in range(arr_u.shape[1]):
    #         velocity_vector[k] = [arr_u[i,j], arr_v[i,j]]
    #         k += 1

    points = np.zeros([len(x)*len(y), 2])
    k = 0

    for i in x:
        for j in y:
            points[k,0] = i
            points[k,1] = j
            k += 1

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
    cells = []

    for i in range(len(tri.simplices)):
        cells.append(Cell(i))

    points_obj = []
    for i in range(len(points)):
        new_point = Point(points[i])
        # new_point.U = velocity_vector[i]
        points_obj.append(new_point)


    # give every cell his neighbors and cells 
    for i in range(len(tri.simplices)):
        c = cells[i]

        bps = [] # List of boundary points for current cell
        for index in tri.simplices[i]:
            bps.append(points_obj[index])
        c.set_boundary_points(bps)
        c.calc_center()
        c.calc_volume()

        for new_neighbor in tri.neighbors[i]:
            if not(new_neighbor == -1):
                c.add_neighbor(cells[new_neighbor])

        # For demonstration give every cell a random value of U 

        # c.assign_random_U()

    if plot_cells:
        plt.triplot(points[:,0], points[:,1], tri.simplices)

        plt.plot(points[:,0], points[:,1], 'o')
        for j, p in enumerate(points):

            plt.text(p[0], p[1], j, ha='right') # label the points

        for j, s in enumerate(tri.simplices):

            p = points[s].mean(axis=0)

            plt.text(p[0], p[1], '#%d' % j, ha='center') # label triangles

        plt.xlim(-0.5, 1.5); plt.ylim(-0.5, 1.5)


    return [cells, points_obj]