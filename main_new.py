# Import Libraries
from numpy import array, true_divide
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

def scatter_quantity(ax, cells):
    x_plot,y_plot,z_plot = get_scatter_values(cells)
    plt.cla()
    ax.scatter(x_plot,y_plot,z_plot)
    ax = plt.gca()
    plt.pause(0.1)

def main():
    
    # numerical parameters
    t = 0
    t_end = 100
    dt = 0.01
    courant_fac = 0.4
    n = 50

    # display parameters
    c1='blue' #blue
    c2='red' #green
    rho_scale = 1.23
    color = False
    scatter = False
    image_paint = True

    # Creating the mesh
    start = timeit.default_timer()
    [cells, points, faces] = create_mesh_rect(n=n, plot_cells=True)
    print(f"Total Cell count {len(cells)}")
    print(f"Mesh Runtime : {timeit.default_timer()-start}")
    # check_cells(cells)

    possible_dts = []
    i=0
    for c in cells:
        if c.number == int(0.5*(n-1)**2):
            c.m, c.mu, c.mv, c.e = getConserved(1.0, 0.01, 0.0, 2.5, c.volume)
        else:
            c.m, c.mu, c.mv, c.e = getConserved(1.0, 0.0, 0.0, 2.5, c.volume)

        c.calc_primitives()

        possible_dts.append(courant_fac * np.min( c.longest_side / \
            (np.sqrt( atm.gamma*c.p/c.rho ) + np.sqrt(c.u**2+c.v**2)) ))
        i+=1
    dt = min(possible_dts)
    print(f"*** Starting timestep dt: {dt}")

    if scatter:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.set_xlabel('x')
        ax.set_ylabel('y')

    i = 0
    while t<t_end:
        

        for c in cells:
            # calculate gradients
            c.calc_gradients_central()

        for c in cells:
            c.extrapol_in_time(dt)
    

        for c in cells:
            c.extrapol2faces()
        tc(cells)
        tf(faces)

        if scatter == True:
            scatter_quantity(ax, cells)
            
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

        
        if color:
            for c in cells:
                u_total = np.sqrt(c.u**2+c.v**2)
                u_total_norm = u_total/0.2
                rho_total_norm = c.rho/rho_scale
                color = colorFader(c1,c2,mix=u_total_norm)
                plt.fill([c.boundary_points[0].X, c.boundary_points[1].X,
                        c.boundary_points[3].X,c.boundary_points[2].X, ], 
                        [c.boundary_points[0].Y, c.boundary_points[1].Y, 
                        c.boundary_points[3].Y, c.boundary_points[2].Y], color)

        


        if image_paint:
            dim = n-1
            A = np.zeros((dim,dim))
            k = 0
            for row in reversed(range(dim)):
                for col in range(dim):
                    u_total = np.sqrt(cells[k].u**2+cells[k].v**2)
                    u_total_norm = u_total/0.2
                    A[row,col] = u_total_norm
                    k +=1
            A = np.true_divide(A,0.2)
        
            plt.cla()
            plt.imshow(A)
            plt.pause(0.5)
            

        # file_name = 'mesh'+f'{i}'+'.pdf'
        # plt.savefig(file_name)

        #subprocess.Popen(['xdg-open '+file_name], shell=True)
        i+=1
    
    plt.show()
        
        


if __name__ == "__main__":
    main()
