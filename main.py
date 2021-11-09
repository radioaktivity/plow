# Import Libraries
from numpy import array
import scipy.spatial as sp
import timeit
import matplotlib.pyplot as plt

from cell import *
from create_mesh import *
from convert import *
from numerical_functions import *
from global_proporties import *


# Import Classes

def get_tangent(x,y,dxdy, width=0.1):
    x_tan = np.linspace(x-width,x+width,10)
    print(x_tan)
    y_tan = dxdy*(x_tan - np.ones([1,10])[0]*x) + np.ones([1, 10])[0]*y
    print(y_tan)
    return [x_tan, y_tan]
    

def main():
    atm.gamma
    t = 0
    t_end = 100
    nx = 3
    ny = 2

    start = timeit.default_timer()
    [cells, points] = create_mesh(nx=nx,ny=ny, plot_cells=True)
    print(f"Mesh Runtime : {timeit.default_timer()-start}")

    xplot = []
    rhoplot = []
    tangents = []
    faces = set()
    for c in cells:
        print(c)
        print('-'*50)
        [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy] = c.calc_gradients_weighted_sum()
        c.extrapol2faces()
        xplot.append(c.center.X)
        rhoplot.append(c.rho)
        tangents.append(get_tangent(c.center.X, c.rho, rho_dx, width=1/(3*nx)))
        for f in c.faces:
            faces.add(f)

    for f in faces:
        plt.scatter(f.center.X, f.rho_L, color=['green'], marker='^')
        plt.scatter(f.center.X, f.rho_R, color=['red'], marker='^')



    print("-"*50)
    plt.scatter(xplot,rhoplot)
    for tan in tangents:
        plt.plot(tan[0], tan[1], 'r')
    plt.show()


if __name__ == "__main__":
    main()
