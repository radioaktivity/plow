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

def text_values_in_cell(cell,gradients=False,primitives=False,fluxes=False):
    text_ofset = atm.linespacing
    if not(atm.printtextincells):
        fluxes =False
        primitives = False
        gradients = False
    if primitives:
        plt.text(cell.center.X, cell.center.Y+text_ofset*4, f'rho:{np.round(cell.rho, decimals=2)}')
        plt.text(cell.center.X, cell.center.Y+text_ofset*3, f'u:{np.round(cell.u, decimals=2)}')
        plt.text(cell.center.X, cell.center.Y+text_ofset*2, f'v:{np.round(cell.v, decimals=2)}')
        plt.text(cell.center.X, cell.center.Y+text_ofset*1, f'p:{np.round(cell.p, decimals=2)}')
    if fluxes:
        plt.text(cell.center.X, cell.center.Y-text_ofset*1, f'm:{np.round(cell.m, decimals=2)}')
        plt.text(cell.center.X, cell.center.Y-text_ofset*2, f'mu:{np.round(cell.mu, decimals=2)}')
        plt.text(cell.center.X, cell.center.Y-text_ofset*3, f'mv:{np.round(cell.mv, decimals=2)}')
        plt.text(cell.center.X, cell.center.Y-text_ofset*4, f'e:{np.round(cell.e, decimals=2)}')
    
    if gradients:
        [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy] = cell.gradients 
        plt.text(cell.center.X, cell.center.Y-text_ofset*1, f'rho_dx:{np.round(rho_dx, decimals=2)}')
        plt.text(cell.center.X, cell.center.Y-text_ofset*2, f'u_dx:{np.round(u_dx, decimals=2)}')
        plt.text(cell.center.X, cell.center.Y-text_ofset*3, f'v_dx:{np.round(v_dx, decimals=2)}')
        plt.text(cell.center.X, cell.center.Y-text_ofset*4, f'p_dx:{np.round(p_dx, decimals=2)}')
        plt.text(cell.center.X, cell.center.Y-text_ofset*5, f'vol:{np.round(cell.volume, decimals=2)}')
        # plt.text(cell.center.X, cell.center.Y-text_ofset*6, f'vec_x:{np.round(cell.ns_neighbor[0][0], decimals=2)}')
        

def text_values_on_face(face, primitives=False, fluxes=False):
    text_ofset = atm.linespacing
    tex_offset_x = 0
    if not(atm.printtextincells):
        fluxes =False
        primitives = False
    if primitives:

        plt.text(face.center.X-tex_offset_x, face.center.Y+text_ofset*4, f'rho_L:{np.round(face.rho_L , decimals=2)}')
        plt.text(face.center.X-tex_offset_x, face.center.Y+text_ofset*3, f'u_L:{np.round(face.u_L , decimals=2)}')
        plt.text(face.center.X-tex_offset_x, face.center.Y+text_ofset*2, f'v_L:{np.round(face.v_L , decimals=2)}')
        plt.text(face.center.X-tex_offset_x, face.center.Y+text_ofset*1, f'p_L:{np.round(face.p_L , decimals=2)}')

        plt.text(face.center.X, face.center.Y, f'Face between {[c.number for c in face.cells_connected]}')

        plt.text(face.center.X+tex_offset_x, face.center.Y+text_ofset*4, f'rho_R:{np.round(face.rho_R , decimals=2)}')
        plt.text(face.center.X+tex_offset_x, face.center.Y+text_ofset*3, f'u_R:{np.round(face.u_R , decimals=2)}')
        plt.text(face.center.X+tex_offset_x, face.center.Y+text_ofset*2, f'v_R:{np.round(face.v_R , decimals=2)}')
        plt.text(face.center.X+tex_offset_x, face.center.Y+text_ofset*1, f'p_R:{np.round(face.p_R , decimals=2)}')

    if fluxes:
        plt.text(face.center.X, face.center.Y+text_ofset*4, f'm X:{np.round(face.flux_Mass , decimals=2)}')
        plt.text(face.center.X, face.center.Y+text_ofset*3, f'mu X:{np.round(face.flux_Momx , decimals=2)}')
        plt.text(face.center.X, face.center.Y+text_ofset*2, f'mv X:{np.round(face.flux_Momy , decimals=2)}')
        plt.text(face.center.X, face.center.Y+text_ofset*1, f'e X:{np.round(face.flux_Energy , decimals=2)}')
    

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
