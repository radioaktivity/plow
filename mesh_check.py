# Import Libraries
from numpy import array
from numpy.core.numeric import convolve
from numpy.matrixlib import defmatrix
import scipy.spatial as sp
import timeit
import matplotlib.pyplot as plt

from cell import *
from create_mesh import *
from convert import *
from numerical_functions import *
from global_proporties import *

def text_values_in_cell(cell):
    text_ofset = 0.012
    plt.text(cell.center.X, cell.center.Y+text_ofset*4, f'rho:{np.round(cell.rho, decimals=2)}')
    plt.text(cell.center.X, cell.center.Y+text_ofset*3, f'u:{np.round(cell.u, decimals=2)}')
    plt.text(cell.center.X, cell.center.Y+text_ofset*2, f'v:{np.round(cell.v, decimals=2)}')
    plt.text(cell.center.X, cell.center.Y+text_ofset*1, f'p:{np.round(cell.p, decimals=2)}')


def check_cells(cells):
    print('-'*50)
    print('Checking Cells...')
    
    bboundary_points = False
    bcenter_exists = False
    bfaces_exist = False
    bneighbors_exist = False
    bvolume = False
    bdis2faces = False
    bns_neighbor = False
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
        if c.ns_neighbor[0][0] == None:
            bns_neighbor = True
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
    elif bns_neighbor:
        raise Exception("check normmal to neighbors")
    elif blongest_side:
        raise Exception("check longest side calculation")
    else:
        print('Cells OK')
