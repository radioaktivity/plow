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

    start = timeit.default_timer()
    [cells, points] = create_mesh(nx=3,ny=2, plot_cells=True)
    print(f"Mesh Runtime : {timeit.default_timer()-start}")

    xplot = []
    rhoplot = []
    tangents = []
    for c in cells:
        print(c)
        print('-'*50)
        [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy] = c.calc_gradients_weighted_sum()
        xplot.append(c.center.X)
        rhoplot.append(c.rho)
        tangents.append(get_tangent(c.center.X, c.rho, rho_dx))

    print("-"*50)
    plt.scatter(xplot,rhoplot)
    for tan in tangents:
        plt.plot(tan[0], tan[1], 'r')
    plt.show()


if __name__ == "__main__":
    main()
