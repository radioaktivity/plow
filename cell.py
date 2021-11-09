import numpy as np
import random

from point import *
from face import *
from convert import *
from global_proporties import *

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
        self.dis2faces = [None,None,None]
        self.dis2neighbors = [None, None, None]
        self.longest_side = 0 # needed for courant number

        # primitives
        self.rho = 0
        self.p = 0
        self.u = 0
        self.v = 0
        self.gradients = None

        # conservatives
        self.m = 0
        self.mu = 0
        self.mv = 0
        self.e = 0


    def calc_primitives(self):
        self.rho, self.u, self.v, self.p = \
            getPrimitive(self.m, self.mu, self.mv, self.e, self.volume)
        

    def calc_all(self):
        if self.boundary_points == []:
            raise Exception('Fatal: No Boundary Points added')
        if self.neighbors == []:
            raise Exception('Fatal: No Neighbors added')
        if self.faces == []:
            raise Exception('Fatal: No faces added')
    
        self.calc_distances_neigbhors()
        self.calc_distances_faces()
        self.calc_volume()


    def set_boundary_points(self, list_of_points:Point):
        self.boundary_points = list_of_points
        self.calc_center()
    
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
            for f in self.faces:
                # Tell the newly created faces, that this cell (self) is connected to it
                f.on_cell(self)
        else: 
            raise Exception("Rectangular cells not yet permitted")

    def face_neighbor_check(self):
        # loops through all neighbors, its faces and compares every one with own faces
        # if it is the same, the neigbors face is taken, if not the own cell is keptcl

        # in all neighbors brows all their faces
        for n in self.neighbors:
            for j in range(len(n.faces)):
                fn = n.faces[j]

                # for one face of the neighbor cell brows all of the own faces
                for i in range(len(self.faces)):
                    f = self.faces[i]
                    equal = f.is_equal_to(fn)

                    # if neighbor cell fn is equal to the own set the own
                    # cell to this cell 
                    if equal:
                        self.faces[i] = fn
                         # tell the face that this cell (self) is connected to it
                        fn.on_cell(self) 

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
        self.longest_side = max([a,b,c])
    
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

            vol.append(n.volume)

            if not(cn[0] == 0): # When 0 centers lie in the same x so no gradient in x-Direction of this cell
                rho_dx += n.volume * (self.rho - n.rho)/(cn[0])
                u_dx += n.volume * (self.u - n.u)/(cn[0])
                v_dx += n.volume * (self.v - n.v)/(cn[0])
                p_dx += n.volume * (self.p - n.p)/(cn[0])

            if not(cn[1]==0):
                rho_dy += n.volume * (self.rho - n.rho)/(cn[1])
                u_dy += n.volume * (self.u - n.u)/(cn[1])
                v_dy += n.volume * (self.v - n.v)/(cn[1])
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

        self.gradients = [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy]

    def extrapol_in_time(self, dt, gamma = 5/3):
        [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy] = self.gradients
        rho_prime = self.rho - 0.5*dt * ( self.u * rho_dx + self.rho * u_dx + self.v * rho_dy + self.rho * v_dy)
        u_prime  = self.u  - 0.5*dt * ( self.u * u_dx + self.v * u_dy + (1/self.rho) * p_dx )
        v_prime  = self.v  - 0.5*dt * ( self.u * v_dx + self.v * v_dy + (1/self.rho) * p_dy )
        p_prime   = self.p   - 0.5*dt * ( gamma*self.p * (u_dx + v_dy)  + self.u * p_dx + self.v * p_dy )

        self.rho = rho_prime
        self.p = u_prime
        self.u = v_prime
        self.v = p_prime


    def extrapol2faces(self):

        [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy] = self.gradients

        for [f,fn] in zip(self.faces, self.dis2faces):
            rho_face = self.rho + rho_dx * fn[0] + rho_dy * fn[1]
            u_face = self.u + u_dx * fn[0] + u_dy * fn[1]
            v_face = self.v + v_dx * fn[0] + v_dy * fn[1]
            p_face = self.p + p_dx * fn[0] + p_dy * fn[1]
            f.get_primitive_value(rho_face, u_face, v_face, p_face)
    
    def get_flux_and_apply(self, dt):

        for f, cnf in zip(self.faces, self.dis2faces):
            '''
            when normal ON face and normal TO face point in one directoin
            substract the flux from this face 
                (same direction -> outside -> negative by definition)
            '''
            direction = np.dot(f.n,cnf)
            direction = - direction/abs(direction) # 1 when flow in, -1 when flow out

            vol = f.get_vol_of_neighbor(self) # get the volume of the current neighbor cell

            self.m += direction * vol * f.m *dt
            self.mu += direction * vol * f.mu * dt
            self.mv += direction * vol * f.mv * dt
            self.e += direction * vol * f.e * dt
        
        # update primitive values
        self.calc_primitives()

    def _str__(self):
        list_neighbors = ''
        for n in self.neighbors:
            list_neighbors = list_neighbors + f';{n.number}\n'
        return f"#Cell {self.number}, #center {self.center.__str__()}\n with the #Neighbors {list_neighbors} " + \
            f"and #boundary points \n {[p.__str__() for p in self.boundary_points]} \n" + \
                f"#Faces: {[f for f in self.faces]}\n"+\
                f"#Faces: {[f.__str__() for f in self.faces]}\n"+\
                f"#volume: {self.volume}\n________________________________________________________________________"
    def __str__(self):
        return f'#Cell {self.number}, '+\
                f'primitives: {self.rho, self.u, self.v, self.p}'+\
                f'gradients: {self.gradients}'+\
                f'conserved: {self.m, self.mu, self.mv, self.e}'+\
                f'\n________________________________________________________________________'

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

    cell1.face_neighbor_check()
    cell2.face_neighbor_check()
    
    cell1.calc_all()
    cell2.calc_all()


    print(cell1.calc_gradients_weighted_sum())

