# Import Libraries
from numpy import array
import scipy.spatial as sp
import timeit
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Import files
from cell import *
from create_mesh import *
from convert import *
from numerical_functions import *
from global_proporties import *
import mesh_check 

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(colors.to_rgb(c1))
    c2=np.array(colors.to_rgb(c2))
    return colors.to_hex((1-mix)*c1 + mix*c2)


def main():
    
    # numerical parameters
    t = 0
    t_end = 100
    dt = 0.01
    nx = 5
    ny = 5
    courant_fac = 0.4

    # display parameters
    c1='blue' #blue
    c2='red' #green
    rho_scale = 1

    # Creating the mesh
    start = timeit.default_timer()
    [cells, points, faces] = create_mesh(nx=nx,ny=ny, plot_cells=True)
    print(f"Total Cell count {len(cells)}")
    print(f"Mesh Runtime : {timeit.default_timer()-start}")
    mesh_check.main(cells, points, faces)

    possible_dts = []
    for c in cells:
        m, mu, mv, e = getConserved(abs(np.random.random(1)[0]), np.random.random(1)[0], 
                                    np.random.random(1)[0], abs(np.random.random(1)[0]),
                                    c.volume)
        c.m = m
        c.mu = mu
        c.mv = mv
        c.e = e
        c.calc_primitives()


        possible_dts.append(courant_fac * np.min( c.longest_side / (np.sqrt( atm.gamma*c.p/c.rho ) + np.sqrt(c.u**2+c.v**2)) ))
    dt = min(possible_dts)


    while t<t_end:

        for c in cells:
            # calculate gradients
            c.calc_gradients_weighted_sum()
            print(c)
            # extrapolate half-step in time

            c.extrapol_in_time(dt)

        # extrapolete in space to face centers

            c.extrapol2faces()

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

        plt.show()
        for c in cells:
            print(f"cell number {c.number} has rho: {c.rho}")
            color = colorFader(c1,c2,mix=c.rho/rho_scale)
            plt.fill([c.boundary_points[0].X, c.boundary_points[1].X, c.boundary_points[2].X], 
                    [c.boundary_points[0].Y, c.boundary_points[1].Y, c.boundary_points[2].Y], color)
            
        plt.show()


if __name__ == "__main__":
    main()
