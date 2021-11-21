import numpy as np
import matplotlib.pyplot as plt

from create_mesh_1d import *
from cell1d import *
from face1d import *
from global_proporties import *

def plot_pressure(cells):
    p = []
    rho = []
    u = []
    pos = []
    for c in cells:
        p.append(c.p)
        u.append(c.u)
        rho.append(c.rho)
        pos.append(c.center.X)
    plt.cla()
    plt.plot(pos,rho, label='rho')
    plt.plot(pos,u, label='u')
    plt.plot(pos,p, label='p')
    plt.legend()
    plt.ylim((-1,3))
    plt.pause(1)

def impulse_initial(cells, n):
    for i, c in enumerate(cells):
        if i==int(n/2):
            c.rho = 1
            c.u = 0.01
            c.p = 2.5
        else:
            c.rho = 1
            c.u = 0
            c.p = 2.5

    return cells

def step_intial(cells, n):
    for i, c in enumerate(cells):
        if (i<=(int(n/2)+int(n/6))) and (i>=(int(n/2)-int(n/6))):
            c.rho = 1
            c.u = 0
            c.p = 2.5
        else:
            c.rho = 1
            c.u = 0.1
            c.p = 2.5

    return cells

def exponential_boundary(cells, n):
    for i, c in enumerate(cells):
        if i<int(n/2):
            c.rho = 1
            c.u = 0.01
            c.p = 2.5
        else:
            c.rho = 1
            c.u = 0
            c.p = 2.5

    return cells

if __name__ == '__main__':

    t = 0
    t_end = 100
    dt = 0
    courant_fac = 0.4
    n = 40

    [cells, faces] = create_mesh(n=n)

    possible_dts = []

    #cells = exponential_boundary(cells)

    for i, c in enumerate(cells):
        c.calc_conserved()
    
    cells = impulse_initial(cells, n)
    plot_pressure(cells)

    for c in cells:
        # Calculate new dt by the courant number in every cell and taking the smallest result
        possible_dts.append( c.distance / (np.sqrt( atm.gamma*c.p/c.rho ) + c.u**2))
    dt = min(possible_dts) * courant_fac
    print(f"dt is {dt}")

    while t<t_end:

        for c in cells:
            c.calc_gradients()

        for c in cells:
            c.extrapol_in_time(dt)

        for c in cells:
            c.extrapol2faces()
        
        for f in faces:
            f.calcFlux()

        for c in cells:
            c.applyFlux(dt)
        

        plot_pressure(cells)

        possible_dts = []
        for c in cells:
            # Calculate new dt by the courant number in every cell and taking the smallest result
            possible_dts.append( c.distance / (np.sqrt( atm.gamma*c.p/c.rho ) + c.u**2) )
        dt = min(possible_dts) * courant_fac
        t += dt
    plt.show()