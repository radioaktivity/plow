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


def main():
    atm.gamma

    start = timeit.default_timer()
    [cells, points] = create_mesh(nx=2,ny=2, plot_cells=True)
    print(f"Mesh Runtime : {timeit.default_timer()-start}")

    for c in cells:
        c.calc_gradients_weighted_sum()

    print("-"*50)

    plt.show()


if __name__ == "__main__":
    main()
