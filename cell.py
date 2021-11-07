import numpy as np

from point import *
from face import *
import gc


class Cell:
    def __init__(self, number):
        # organizing
        self.number = number
        self.neighbors = []
        self.faces = []

        # geometry
        self.center = Point()
        self.boundary_points = []
        self.volume = 0
        self.dis2faces = []
        self.dis2neighbors = []

        # primitive values (at the center)
        self.rho = 0
        self.p = 0
        self.u = 0
        self.v = 0



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
        for n in self.neighbors:
            self.center.distance(n.center)

    def calc_distances_faces(self):
        # Calculates distances to all face centers of the cell
        for f in self.faces:
            self.center.distance(f.center)

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

