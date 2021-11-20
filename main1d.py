import numpy as np
import matplotlib.pyplot as plt

from create_mesh_1d import *
from cell1d import *
from face1d import *


if __name__ == '__main__':

    t = 0
    t_end = 100
    dt = 0
    courant_fac = 0.4
    n = 50

    [cells, faces] = create_mesh(n=n)

    possible_dts = []

    for c in cells:
        # Calculate new dt by the courant number in every cell and taking the smallest result
        possible_dts.append( c.distance / (np.sqrt( atm.gamma*c.p/c.rho ) + c.u**2))
    dt = min(possible_dts) * courant_fac


    while t<t_end:

        for c in cells:
            c.calc_gradients()

        for c in cells:
            c.extrapol_in_time(dt)
        
        for c in cells:
            c.extrapol2faces()
        
        for f in faces:
            f.getFlux()

        for c in cells:
            c.applyFlux(dt)
        

        possible_dts = []
        for c in cells:
            # Calculate new dt by the courant number in every cell and taking the smallest result
            possible_dts.append( c.distance / (np.sqrt( atm.gamma*c.p/c.rho ) + c.u**2) )
        dt = min(possible_dts) * courant_fac
        t += dt