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
    [cells, points] = create_mesh(nx=300,ny=300, plot_cells=False)
    print(f"Mesh Runtime : {timeit.default_timer()-start}")
    print("-"*50)
    


if __name__ == "__main__":
    main()
