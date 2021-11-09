# Import Libraries
from numpy import array
import scipy.spatial as sp
import timeit
import matplotlib.pyplot as plt

# Import files
from cell import *
from create_mesh import *
from convert import *
from numerical_functions import *
from global_proporties import *
import mesh_check 



def main():
    atm.gamma
    t = 0
    t_end = 100
    dt = 0.01
    nx = 3
    ny = 2
    courant_fac = 0.4

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
            print(f)

        possible_dts = []

        # apply fluxes 
        for c in cells:
            c.get_flux_and_apply(dt)      

        # Calculate new dt by the courant number in every cell and taking the smallest result
            possible_dts.append(courant_fac * np.min( c.longest_side / (np.sqrt( atm.gamma*c.p/c.rho ) + np.sqrt(c.u**2+c.v**2)) ))
        dt = min(possible_dts)

        # update time
        t += dt
    



if __name__ == "__main__":
    main()
