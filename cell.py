import numpy as np

from point import *
from face import *


class Cell:
    def __init__(self, number):
        # organizing
        self.number = number
        self.neighbors = []

        # geometry
        self.center = 0
        self.boundary_points = []
        self.volume = 0

        # primitive values
        self.rho = 0
        self.p = 0
        self.u = 0
        self.v = 0


    def set_boundary_points(self, list_of_points:Point):
        self.boundary_points = list_of_points
    
    def add_neighbor(self, new_neighbors):
        self.neighbors.append(new_neighbors)

    def calc_center(self):
        boundary_points_ko =np.zeros([len(self.boundary_points),2])
        k = 0
        for p in self.boundary_points:
            boundary_points_ko[k] = p.getValue()
            k+=1
        self.center  = boundary_points_ko.mean(axis=0)

    # def assign_random_U(self):
    #     u_from_points = []
    #     v_from_points = []
    #     for p in self.boundary_points:
    #         u_from_points.append(p.U[0])
    #         v_from_points.append(p.U[1])
    #     self.U = [np.array(u_from_points).mean(),np.array(v_from_points).mean()]

    def calc_volume(self):
        a = self.boundary_points[0].distance(self.boundary_points[1])
        b = self.boundary_points[0].distance(self.boundary_points[2])
        surface = a * b * 0.5 
        self.volume = surface  

    def __str__(self):
        list_neighbors = ''
        for n in self.neighbors:
            list_neighbors = list_neighbors + f';{n.number}\n'
        return f"Cell at {self.center}\n with the neighbors {list_neighbors} " + \
            f"and boundary points \n {[p.__str__() for p in self.boundary_points]}  \n value U {self.U}\n" + \
                f"volume: {self.volume}\n__________________"

if __name__ == "__main__":
    cell1 = Cell(0)
    print(cell1.center) 