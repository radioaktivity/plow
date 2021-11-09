# Import Libraries
from numpy import array
from numpy.core.numeric import convolve
import scipy.spatial as sp
import timeit
import matplotlib.pyplot as plt

from cell import *
from create_mesh import *
from convert import *
from numerical_functions import *
from global_proporties import *


def get_tangent(x,y,dxdy, width=0.1):
    x_tan = np.linspace(x-width,x+width,10)
    print(x_tan)
    y_tan = dxdy*(x_tan - np.ones([1,10])[0]*x) + np.ones([1, 10])[0]*y
    print(y_tan)
    return [x_tan, y_tan]

def plot_tangents(cells):
    xplot = []
    rhoplot = []
    tangents = []
    faces = set()

    for c in cells:
        print(c)
        print('-'*50)
        [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy] = c.calc_gradients_weighted_sum()
        c.extrapol2faces()
        xplot.append(c.center.X)
        rhoplot.append(c.rho)
        tangents.append(get_tangent(c.center.X, c.rho, rho_dx, width=1/2 * c.longest_side))
        for f in c.faces:
            faces.add(f)

    for f in faces:
        plt.scatter(f.center.X, f.rho_L, color=['green'], marker='^')
        plt.scatter(f.center.X, f.rho_R, color=['red'], marker='^')



    print("-"*50)
    plt.scatter(xplot,rhoplot)
    for tan in tangents:
        plt.plot(tan[0], tan[1], 'r')
    plt.show()


def check_cells(cells):
    print('-'*50)
    print('Checking Cells...')
    
    bboundary_points = False
    bcenter_exists = False
    bfaces_exist = False
    bneighbors_exist = False
    bvolume = False
    bdis2faces = False
    bdis2neighbors = False
    blongest_side = False

    for c in cells:
        if c.boundary_points == []:
            bboundary_points = True
        if not(c.center.isDefined()):
            bcenter_exists = True
        if c.faces == []:
            bfaces_exist = True
        if c.neighbors == []:
            bneighbors_exist = True
        if c.volume == 0:
            bvolume = True
        if c.dis2neighbors[0][0] == None:
            bdis2neighbors = True
        if c.dis2faces[0][0] == None:
            bdis2faces = True
        if c.longest_side == 0:
            blongest_side = True

    if bcenter_exists:
        raise Exception("check centers")
    elif bfaces_exist:
        raise Exception("check faces")
    elif bneighbors_exist:
        raise Exception("check neighbors")
    elif bboundary_points:
        raise Exception("check boundary points")
    elif  bvolume:
        raise Exception("check volume")
    elif bdis2faces:
        raise Exception("check normal to faces")
    elif bdis2neighbors:
        raise Exception("check normmal to neighbors")
    elif blongest_side:
        raise Exception("check longest side calculation")
    else:
        print('Cells OK')


def main(cells, points, faces):
    check_cells(cells)