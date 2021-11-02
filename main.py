# Import Libraries
from numpy import array
import scipy.spatial as sp
from cell import *
from create_mesh import *
from convert import *
from numerical_functions import *
from global_proporties import *
import matplotlib.pyplot as plt

# Import Classes


def main():
    atm.gamma

    [cells, points] = create_mesh(nx=2,ny=2, plot_cells=True)
    for c in cells:
        print(c)
    x = np.array([cell.center[0] for cell in cells])
    y = np.array([cell.center[1] for cell in cells])
    
    u = np.array([cell.U[0] for cell in cells])
    v = np.array([cell.U[1] for cell in cells])

    print(x.shape,y.shape,u.shape,v.shape)


    u_abs = np.sqrt(np.multiply(u,u)*np.multiply(v,v))

    plt.scatter(x,y,c=u_abs)
    plt.show()
    


if __name__ == "__main__":
    main()
