import numpy as np
import matplotlib.pyplot as plt

from create_mesh_1d import *
from cell1d import *
from face1d import *
from global_proporties import *


def plot_pressure(cells, title=' ', pause=1, faces=None):
    p = []
    rho = []
    u = []
    pos = []

    dim = cells[0].volume
    for c in cells:
        p.append(c.p)
        u.append(c.u)
        rho.append(c.rho)
        pos.append(c.center.X)
    plt.cla()
    plt.plot(pos,rho, label='rho')
    plt.plot(pos,u, label='u')
    plt.plot(pos,p, label='p')


    if not(faces==None):
        u_L = []
        u_R = []
        pos_F = []
        for f in faces:
            u_L.append(f.u_L)
            u_R.append(f.u_R)
            pos_F.append(f.center.X)

        plt.scatter(pos_F-np.ones(len(pos_F))*dim/5, u_L, label='u_L')
        plt.scatter(pos_F+np.ones(len(pos_F))*dim/5, u_R, label='u_R')


    plt.legend()
    plt.title(title)
    plt.ylim((-1,4))
    plt.pause(pause)

def impulse_initial(cells, n, size=0.3):
    for i, c in enumerate(cells):
        if i==int(n/2):
            c.rho = 1
            c.u = size
            c.p = 2.5
        else:
            c.rho = 1
            c.u = 0
            c.p = 2.5

    return cells

def step_intial(cells, n, size=0.1):
    for i, c in enumerate(cells):
        if (i<=(int(n/2)+int(n/6))) and (i>=(int(n/2)-int(n/6))):
            c.rho = 1
            c.u = 0
            c.p = 2.5
        else:
            c.rho = 1
            c.u = size
            c.p = 2.5

    return cells

def exponential_boundary(cells, n):
    u0 =1 
    for i, c in enumerate(cells):
        if i<int(n/2):
            c.rho = 1
            c.u = u0 * i/n
            c.p = 2.5
        else:
            c.rho = 1
            c.u = u0 * (1-i/n)
            c.p = 2.5

    return cells

if __name__ == '__main__':

    t = 0
    t_end = 10
    dt = 0
    pause = 0.01
    nth_plot = 100
    courant_fac = 0.4
    n = 50

    [cells, faces] = create_mesh(n=n)

    possible_dts = []

    #cells = exponential_boundary(cells)
    
    # cells = exponential_boundary(cells, n)
    cells = step_intial(cells, n, size=-1.5)
    cells = impulse_initial(cells, n, size=3)
    plot_pressure(cells)


    for i, c in enumerate(cells):
        c.calc_conserved()

    for c in cells:
        # Calculate new dt by the courant number in every cell and taking the smallest result
        possible_dts.append( c.distance / (np.sqrt( atm.gamma*c.p/c.rho ) + c.u**2))
    dt = min(possible_dts) * courant_fac
    print(f"dt is {dt}")

    count = 0
    while t<t_end:
        if count%nth_plot==0:
            print(t)
            plot_pressure(cells, title='flux applied', pause=pause, faces=None)

        for c in cells:
            c.calc_gradients(type='central')

        for c in cells:
            c.extrapol_in_time(dt)
        #plot_pressure(cells, title='in time')
        for c in cells:
            c.extrapol2faces()
        #plot_pressure(cells, title='to face')
        for f in faces:
            f.calcFlux()
        
        for c in cells:
            c.applyFlux(dt)
        

    

        possible_dts = []
        for c in cells:
            # Calculate new dt by the courant number in every cell and taking the smallest result
            possible_dts.append( c.distance / (np.sqrt( atm.gamma*c.p/c.rho ) + c.u**2) )
        dt = min(possible_dts) * courant_fac
        t += dt
        count += 1
    plt.show()