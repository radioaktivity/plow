# Import Libraries
from numpy import array
import scipy.spatial as sp
import timeit
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import subprocess


# Import files
from cell import *
from create_mesh import *
from convert import *
from numerical_functions import *
from global_proporties import *
from mesh_check import *
from vector_alg import *

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    if mix>1:
        mix = 1
    if mix<0:
        mix = 0
    c1=np.array(colors.to_rgb(c1))
    c2=np.array(colors.to_rgb(c2))
    return colors.to_hex((1-mix)*c1 + mix*c2)

def print_cell(cell, number=10):
    if (cell.number == number):
        print('_'*50)
        print("Cell 10 Values:")
        print(cell.u)
        print(cell.v)
        print('Neighbors')
        for n in cell.neighbors:
            print('Cell', n.number)
            print(n.u)
            print(n.v)
            print('-'*30)


def main():
    
    # numerical parameters
    t = 0
    t_end = 100
    dt = 0.01
    nx = 3
    ny = 3
    courant_fac = 0.4

    # display parameters
    c1='blue' #blue
    c2='red' #green
    rho_scale = 1.23

    # Creating the mesh
    start = timeit.default_timer()
    [cells, points, faces] = create_mesh(nx=nx,ny=ny, plot_cells=True)
    print(f"Total Cell count {len(cells)}")
    print(f"Mesh Runtime : {timeit.default_timer()-start}")
    check_cells(cells)
   
       

    possible_dts = []
    i=0
    for c in cells:
        if c.number in [0,1,2,3,4,5,6,7]:
            c.m, c.mu, c.mv, c.e = getConserved(1.225, 5, 2, 100, c.volume)
        else:
            c.m, c.mu, c.mv, c.e = getConserved(1.225, 5, 0, 100, c.volume)

        c.calc_primitives()

        possible_dts.append(courant_fac * np.min( c.longest_side / (np.sqrt( atm.gamma*c.p/c.rho ) + np.sqrt(c.u**2+c.v**2)) ))
        i+=1
    dt = min(possible_dts)
    print(f"*** Starting timestep dt: {dt}")


    while t<t_end:

        for c in cells:
            # calculate gradients
            c.calc_gradients_weighted_sum()

        for c in cells:
            c.extrapol_in_time(dt)
        

        for c in cells:
            c.extrapol2faces()

        # compute fluxes

        for f in faces:
            f.getFlux()

        possible_dts = []


        for c in cells:
            text_values_in_cell(c)

        for f in faces:
            text_values_on_face(f)

        
        # apply fluxes 
        for c in cells:
            c.get_flux_and_apply(dt)      

        # Calculate new dt by the courant number in every cell and taking the smallest result
            possible_dts.append(courant_fac * np.min( c.longest_side / (np.sqrt( atm.gamma*c.p/c.rho ) + np.sqrt(c.u**2+c.v**2)) ))
        dt = min(possible_dts)

        # update time
        t += dt

        for c in cells:
            print(f"cell number {c.number} has rho: {c.rho}")
            color = colorFader(c1,c2,mix=c.rho/rho_scale)
            plt.fill([c.boundary_points[0].X, c.boundary_points[1].X, c.boundary_points[2].X], 
                    [c.boundary_points[0].Y, c.boundary_points[1].Y, c.boundary_points[2].Y], color)
        
        
        file_name = 'mesh'+f'{np.round(t,decimals=3)}'+'.pdf'
        plt.savefig(file_name)

        # subprocess.Popen(['xdg-open '+file_name], shell=True)
        
    
        


if __name__ == "__main__":
    main()
