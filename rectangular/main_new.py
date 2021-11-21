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

def plot_image(cells, n, pause=0.1):
    dim = n-1
    A = np.zeros((dim,dim))
    k = 0
    for row in reversed(range(dim)):
        for col in range(dim):
            u_total = np.sqrt(cells[k].u**2+cells[k].v**2)
            u_total_norm = u_total/0.2
            A[row,col] = cells[k].p
            k +=1
    # A = np.true_divide(A,0.2)

    plt.cla()
    plt.imshow(A)
    plt.pause(pause)

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



def get_scatter_values(cells):
    x_plot = []
    y_plot = []
    z_plot = []
    for c in cells:

        x_plot.append(c.center.X)
        y_plot.append(c.center.Y)
        z_plot.append(c.p)

    return x_plot, y_plot, z_plot

def impulse_initial(cells, n, size=(0.1, 0)):
    for i, c in enumerate(cells):
        if i==int(n**2/2):
            c.rho = 1
            c.u = size[0]
            c.v = size[1]
            c.p = 2.5
        else:
            c.rho = 1
            c.u = 0
            c.v = 0
            c.p = 2.5

    return cells


def step_intial(cells, n, size=0.1):
    for i, c in enumerate(cells):
        if (i<=(int(n/2)+int(n/6))) and (i>=(int(n/2)-int(n/6))):
            c.rho = 1
            c.u = 0
            c.v = 0
            c.p = 2.5
        else:
            c.rho = 1
            c.u = size
            c.v = 0
            c.p = 2.5
    return cells

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
    courant_fac = 0.2
    n = 16
    nth_turn = 5
    pause = 0.1
    scatter = False
    image_paint = True

    # Creating the mesh
    start = timeit.default_timer()

    [cells, points, faces] = create_mesh_rect(n=n, plot_cells=False)
    print(f"Total Cell count {len(cells)}")
    print(f"Mesh Runtime : {timeit.default_timer()-start}")
    # check_cells(cells)


    # cells = exponential_boundary(cells, n)
    cells = impulse_initial(cells, n, size=(0,0.1))

    for i, c in enumerate(cells):
        c.calc_conserved()

    possible_dts = []
    for c in cells:
        possible_dts.append(np.min( c.longest_side / \
            (np.sqrt( atm.gamma*c.p/c.rho ) + np.sqrt(c.u**2+c.v**2)) ))
    dt = min(possible_dts) * courant_fac
    print(f"*** Starting timestep dt: {dt}")

    plot_image(cells, n, pause=1)
    if scatter:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.set_xlabel('x')
        ax.set_ylabel('y')

    i = 0
    while t<t_end:

        for c in cells:
            c.calc_gradients_central()

        for c in cells:
            c.extrapol_in_time(dt)
            
        for c in cells:
            c.extrapol2faces()

        for f in faces:
            f.getFlux()
            

        # apply fluxes 
        for c in cells:
            c.get_flux_and_apply(dt)      

            if c.rho <= 0:
                raise Exception('Negative Rho')
            if c.p <= 0:
                raise Exception('Negative P')

        possible_dts = []
        for c in cells:
        # Calculate new dt by the courant number in every cell and taking the smallest result
            possible_dts.append(( c.longest_side / (np.sqrt( atm.gamma*c.p/c.rho ) + np.sqrt(c.u**2+c.v**2)) ))
        dt = min(possible_dts) * courant_fac

                        
        if image_paint and (i%nth_turn==0):
            print(t)
            plot_image(cells, n, pause=pause)
        elif scatter and (i%nth_turn==0):
            scatter_quantity(ax, cells)


        # update time
        t += dt

        i+=1
    
    plt.show()
        
        


if __name__ == "__main__":
    main()
