import numpy as np
import random

from point import *
from face import *
from convert import *


class Cell:
    def __init__(self, number):

        # organizing
        self.number = number
        self.neighbors = []
        self.faces = []

        # geometry
        self.center = Point()
        self.boundary_points = []
        self.volume = 0 # or surface in 2D
        self.dis2faces = [np.zeros([2,1]),np.zeros([2,1]),np.zeros([2,1])]
        self.dis2neighbors = [np.zeros([2,1]), np.zeros([2,1]), np.zeros([2,1])]

        # primitive values (at the center)
        self.rho = random.uniform(0,1)
        self.p = random.uniform(0,1)
        self.u = random.uniform(0,1)
        self.v = random.uniform(0,1)



    def set_boundary_points(self, list_of_points:Point):
        self.boundary_points = list_of_points
    
    def add_neighbor(self, new_neighbors):
        self.neighbors.append(new_neighbors)
    
    def create_faces(self):
        # creating all faces out of boundary points

        if len(self.boundary_points) == 3:
            p1 = self.boundary_points[0]
            p2 = self.boundary_points[1]
            p3 = self.boundary_points[2]
            self.faces.append(Face(p1,p2))
            self.faces.append(Face(p3,p2))
            self.faces.append(Face(p1,p3))
        else: 
            raise Exception("Rectangular cells not yet permitted")

    def face_neighbor_check(self):
        # loops through all neighbors, its faces and compares every one with own faces
        # if it is the same, the neigbors face is taken, if not the own cell is keptcl

        for n in self.neighbors:

            for j in range(len(n.faces)):
                fn = n.faces[j]
                for i in range(len(self.faces)):
                    f = self.faces[i]
                    equal = f.is_equal_to(fn)
                    if equal:
                        self.faces[i] = fn

    def calc_center(self):
        # Calculating the center by the mean of all boundary points
        boundary_points_ko =np.zeros([len(self.boundary_points),2])
        k = 0
        for p in self.boundary_points:
            boundary_points_ko[k] = p.getValue()
            k+=1
        center = boundary_points_ko.mean(axis=0)
        self.center = Point(center)

    def calc_volume(self):
        # Calculating the Volume (3D)/ the Surface (2D) of the cell
        a = self.boundary_points[0].distance(self.boundary_points[1])
        b = self.boundary_points[0].distance(self.boundary_points[2])
        c = self.boundary_points[1].distance(self.boundary_points[2])
        s = 0.5*(a+b+c)
        surface = np.sqrt(s*(s-a)*(s-b)*(s-c))
        self.volume = surface  
    
    def calc_distances_neigbhors(self):
        # Calculates distances to all neigbors of the cell
        i = 0
        for n in self.neighbors:
            self.dis2neighbors[i] = self.center.getVecBetween(n.center)
            i+=1

    def calc_distances_faces(self):
        # Calculates distances to all face centers of the cell
        i = 0
        for f in self.faces:
            self.dis2faces[i] = self.center.getVecBetween(f.center)
            i+=1

    def calc_gradients_weighted_sum(self):
        '''
        gradients are calculated by taking the gradient
        from this cell to every neighbor cell and building 
        an avarage which is weighted by the cell volumes
        vol         is Array of neighbor cell volums
        n           is the current neighbor cell
        cn          is the vector between the current neighbor cell and this cell
        drho_dx     is the gradient of rho in the x direction
        '''
        vol = []
        
        [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy]=\
            [0.,0.,0.,0.,0.,0.,0.,0.]

        for [n,cn] in zip(self.neighbors, self.dis2neighbors):
            print(cn)
            vol.append(n.volume)

            rho_dx += n.volume * (self.rho - n.rho)/(cn[0])
            rho_dy += n.volume * (self.rho - n.rho)/(cn[1])

            u_dx += n.volume * (self.u - n.u)/(cn[0])
            u_dy += n.volume * (self.u - n.u)/(cn[1])

            v_dx += n.volume * (self.v - n.v)/(cn[0])
            v_dy += n.volume * (self.v - n.v)/(cn[1])

            p_dx += n.volume * (self.p - n.p)/(cn[0])
            p_dy += n.volume * (self.p - n.p)/(cn[1])


        v_tot = np.sum(vol)

        rho_dx *= 1/v_tot
        rho_dy *= 1/v_tot

        u_dx *= 1/v_tot
        u_dy *= 1/v_tot

        v_dx *= 1/v_tot
        v_dy *= 1/v_tot

        p_dx *= 1/v_tot
        p_dy *= 1/v_tot
        return [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy]

    def __str__(self):
        list_neighbors = ''
        for n in self.neighbors:
            list_neighbors = list_neighbors + f';{n.number}\n'
        return f"#Cell {self.number}, #center {self.center.__str__()}\n with the #Neighbors {list_neighbors} " + \
            f"and #boundary points \n {[p.__str__() for p in self.boundary_points]} \n" + \
                f"#Faces: {[f for f in self.faces]}\n"+\
                f"#Faces: {[f.__str__() for f in self.faces]}\n"+\
                f"#volume: {self.volume}\n__________________"

if __name__ == "__main__":

    cell1 = Cell(0)
    cell2 = Cell(1)

    p1 = Point(np.array([0,0]))
    p2 = Point(np.array([0,1]))
    p3 = Point(np.array([1,1]))
    p4 = Point(np.array([1,0]))

    cell1.set_boundary_points([p1,p2,p3])
    cell2.set_boundary_points([p1,p3,p4])

    cell1.add_neighbor(cell2)
    cell2.add_neighbor(cell1)

    cell1.create_faces()
    cell2.create_faces()
    print("################# BEFORE ##################")
    print(cell1)
    print(cell2)

    cell1.face_neighbor_check()
    cell2.face_neighbor_check()
    print("################# AFTER ##################")
    print(cell1)
    print(cell2)

