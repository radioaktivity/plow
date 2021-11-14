# Import Libraries
from numpy import array
import scipy.spatial as sp
import timeit
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import subprocess
from mpl_toolkits.mplot3d import Axes3D


# Import files
from cell_new import *
from create_mesh_rectangular import *
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

def tc(cells):
    for c in cells:
        text_values_in_cell(c,primitives=True, gradients=True)
def tf(faces):
    for f in faces:
        text_values_on_face(f, primitives=True)

def get_scatter_values(cells):
    x_plot = []
    y_plot = []
    z_plot = []
    for c in cells:

        x_plot.append(c.center.X)
        y_plot.append(c.center.Y)
        z_plot.append(c.u)

    return x_plot, y_plot, z_plot

def main():
    
    # numerical parameters
    t = 0
    t_end = 100
    dt = 0.01
    courant_fac = 0.2
    n = 5

    # display parameters
    c1='blue' #blue
    c2='red' #green
    rho_scale = 1.23
    color = False
    scatter = True

    # Creating the mesh
    start = timeit.default_timer()
    [cells, points, faces] = create_mesh_rect(n=n, plot_cells=True)
    print(f"Total Cell count {len(cells)}")
    print(f"Mesh Runtime : {timeit.default_timer()-start}")
    # check_cells(cells)

    possible_dts = []
    i=0
    for c in cells:
        if c.number in range(0,int((n-1)**2*0.5)):
            c.m, c.mu, c.mv, c.e = getConserved(1, 0.1, 0.0, 2.5, c.volume)
        else:
            c.m, c.mu, c.mv, c.e = getConserved(1, 0, 0, 2.5, c.volume)

        c.calc_primitives()

        possible_dts.append(courant_fac * np.min( c.longest_side / (np.sqrt( atm.gamma*c.p/c.rho ) + np.sqrt(c.u**2+c.v**2)) ))
        i+=1
    dt = min(possible_dts)
    print(f"*** Starting timestep dt: {dt}")
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    i = 0
    while t<t_end:
        
        x_plot,y_plot,z_plot = get_scatter_values(cells)
        plt.cla()
        ax.scatter(x_plot,y_plot,z_plot)
        ax = plt.gca()
        plt.pause(1)

        for c in cells:
            # calculate gradients
            c.calc_gradients_upwind()

        for c in cells:
            c.extrapol_in_time(dt)
        
        x_plot,y_plot,z_plot = get_scatter_values(cells)
        plt.cla()
        ax.scatter(x_plot,y_plot,z_plot)
        ax = plt.gca()
        plt.pause(1)

        for c in cells:
            c.extrapol2faces()
        tc(cells)
        tf(faces)

        # compute fluxes

        for f in faces:
            f.getFlux()


        possible_dts = []

        # apply fluxes 
        for c in cells:
            c.get_flux_and_apply(dt)      

        

        # Calculate new dt by the courant number in every cell and taking the smallest result
            possible_dts.append(courant_fac * np.min( c.longest_side / (np.sqrt( atm.gamma*c.p/c.rho ) + np.sqrt(c.u**2+c.v**2)) ))
        dt = min(possible_dts)

        # update time
        t += dt

        

        for c in cells:
            # print(f"cell number {c.number} has rho: {c.rho} has u: {c.u}")

            if color:
                u_total = np.sqrt(c.u**2+c.v**2)
                u_total_norm = u_total/0.2
                rho_total_norm = c.rho/rho_scale
                color = colorFader(c1,c2,mix=u_total_norm)
                plt.fill([c.boundary_points[0].X, c.boundary_points[1].X,
                        c.boundary_points[3].X,c.boundary_points[2].X, ], 
                        [c.boundary_points[0].Y, c.boundary_points[1].Y, 
                        c.boundary_points[3].Y, c.boundary_points[2].Y], color)

        
        # file_name = 'mesh'+f'{i}'+'.pdf'
        # plt.savefig(file_name)

        #plt.show()
        #ubprocess.Popen(['xdg-open '+file_name], shell=True)
        i+=1
        
        x_plot,y_plot,z_plot = get_scatter_values(cells)
        plt.cla()
        ax.scatter(x_plot,y_plot,z_plot)
        ax = plt.gca()
        plt.pause(1)
        
    
    plt.show()
        
        


if __name__ == "__main__":
    main()
