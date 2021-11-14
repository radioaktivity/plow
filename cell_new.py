import numpy as np
import random

from numpy import linalg

from point import *
from face_new import *
from convert import *
from global_proporties import *
from vector_alg import *

class Cell:
    def __init__(self, number):

        # organizing
        self.number = number
        self.neighbors = []
        self.faces = []
        self.sides = []

        # geometry
        self.center = Point()
        self.boundary_points = []
        self.volume = 0 # or surface in 2D
        self.vecs2faces = []
        self.vecs2neighbors = []
        self.longest_side = 0 # needed for courant number
        self.non_orto_angles = [] # goes along with neighbors

        # primitives
        self.rho = 0
        self.u = 0
        self.v = 0
        self.p = 0
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

        self.calc_volume()
        self.calc_vecs2faces()

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
            p1 = self.boundary_points[0]
            p2 = self.boundary_points[1]
            p3 = self.boundary_points[2]
            p4 = self.boundary_points[3]

            # NOTE : Faces are sorted as: Left, Right, Top, Bottom
            self.faces.append(Face(p3,p1,type='X'))
            self.faces.append(Face(p2,p4,type='X'))
            self.faces.append(Face(p4,p3,type='Y'))
            self.faces.append(Face(p1,p2,type='Y'))
            
            
            for f in self.faces:
                # Tell the newly created faces, that this cell (self) is connected to it
                f.on_cell(self)

    def face_neighbor_check(self):
        '''
        # loops through all neighbors, its faces and compares every one with own faces
        # if it is the same, the neigbors face is taken, if not the own cell is keptcl
        n       is the current neighbor
        fn      is the current neighbor face (index j)
        f       is the current face of cell (index i)
        '''

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

    def face_neighbor_check2(self):
        '''
        Only for tetra-mesh
        takes adwantage of faces being sorted (left right top bottom)
        and neighbors being sorted (left right top bottom)

        ATTENTION : MIXES UP KOORDINATES
        '''
        self.faces[0] = self.neighbors[0].faces[1]
        self.faces[1] = self.neighbors[1].faces[0]
        self.faces[2] = self.neighbors[2].faces[3]
        self.faces[3] = self.neighbors[3].faces[2]
        
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
        if(len(self.boundary_points)):
            # Calculating the Volume (3D)/ the Surface (2D) of the cell
            a = self.boundary_points[0].distance(self.boundary_points[1])
            b = self.boundary_points[0].distance(self.boundary_points[2])
            c = self.boundary_points[1].distance(self.boundary_points[2])
            s = 0.5*(a+b+c)
            surface = np.sqrt(s*(s-a)*(s-b)*(s-c))
            self.volume = surface  
            self.longest_side = max([a,b,c])
        else:
            surface = 0
            # Calculating the Volume (3D)/ the Surface (2D) of the cell
            a = self.faces[0].surface
            b = self.faces[2].surface
            c = np.linalg.norm(self.boundary_points[2].getVecBetween(self.boundary_points[3]))
            s = 0.5*(a+b+c)
            surface += np.sqrt(s*(s-a)*(s-b)*(s-c))
            # Calculating the Volume (3D)/ the Surface (2D) of the cell
            a = self.faces[1].surface
            b = self.faces[3].surface
            c = c
            s = 0.5*(a+b+c)
            surface += np.sqrt(s*(s-a)*(s-b)*(s-c))

            self.longest_side = max([f.surface for f in self.faces])
            self.volume = surface  

    def calc_vecs2faces(self):
        # Calculates distances to all face centers of the cell
        
        for f in self.faces:
            self.vecs2faces.append(self.center.getVecBetween(f.center))
    

    def calc_gradients_upwind(self):
       
        [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy]=\
            [0.,0.,0.,0.,0.,0.,0.,0.]

        # X Direction Gradient 
        if self.u > 0: # Left neighbor is needed
           n = self.neighbors[0] 
           direction  = -1
        else: # Right neighbor is needed
           n = self.neighbors[1]
           direction  = 1

        d12 = self.center.getVecBetween(n.center)
        d = d12[0] # vector to neihbor

        # top face length 
        d = self.faces[2].surface *direction

        rho_dx += (n.rho-self.rho)/d
        u_dx +=  (n.u-self.u)/d
        v_dx +=  (n.v-self.v)/d
        p_dx +=  (n.p-self.p)/d
            

        # Y Direction Gradient 
        if self.v > 0: # bottom neighbor is needed
            n = self.neighbors[3]
            direction = -1
        else: # top neighbor is needed
            n = self.neighbors[2]
            direction = 1

        d12 = self.center.getVecBetween(n.center)
        d = d12[1] # vector to neighbor
        # left face length
        d = self.faces[0].surface *direction

        rho_dy += (n.rho-self.rho)/d
        u_dy +=  (n.u-self.u)/d
        v_dy +=  (n.v-self.v)/d
        p_dy += (n.p-self.p)/d
            
        self.gradients = [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy]


    def extrapol_in_time(self, dt):
        gamma = atm.gamma
        [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy] = self.gradients
        rho, u, v, p = self.rho, self.u, self.v, self.p

        # extrapolate half-step in time
        rho_prime = rho - 0.5*dt * ( u * rho_dx + rho * u_dx + v * rho_dy + rho * v_dy)
        u_prime = u - 0.5*dt * ( u * u_dx + v * u_dy + (1/rho) * p_dx )
        v_prime = v - 0.5*dt * ( u * v_dx + v * v_dy + (1/rho) * p_dy )
        p_prime = p - 0.5*dt * ( gamma*p * (u_dx + v_dy)  + u * p_dx + v * p_dy )
        
        self.rho = rho_prime
        self.u = u_prime
        self.v = v_prime
        self.p = p_prime

    def extrapol2faces(self):

        [rho_dx, rho_dy, u_dx, u_dy, v_dx, v_dy, p_dx, p_dy] = self.gradients


        #left right top bottom
        # x x y y 
        # LEFT
        f = self.faces[0]

        fn = self.center.getVecBetween(f.center)
        
        rho_face = self.rho + rho_dx * fn[0] 
        u_face = self.u + u_dx * fn[0] 
        v_face = self.v + v_dx * fn[0] 
        p_face = self.p + p_dx * fn[0] 

        f.get_primitive_values(rho_face, u_face, v_face, p_face)

        # RIGHT
        f = self.faces[1]

        fn = self.center.getVecBetween(f.center)
        
        rho_face = self.rho + rho_dx * fn[0] 
        u_face = self.u + u_dx * fn[0] 
        v_face = self.v + v_dx * fn[0] 
        p_face = self.p + p_dx * fn[0] 

        f.get_primitive_values(rho_face, u_face, v_face, p_face)

        # TOP
        f = self.faces[2]

        fn = self.center.getVecBetween(f.center)
        
        rho_face = self.rho + rho_dy * fn[1] 
        u_face = self.u + u_dy * fn[1] 
        v_face = self.v + v_dy * fn[1] 
        p_face = self.p + p_dy * fn[1] 

        f.get_primitive_values(rho_face, u_face, v_face, p_face)


        # TOP
        f = self.faces[3]

        fn = self.center.getVecBetween(f.center)
        
        rho_face = self.rho + rho_dy * fn[1] 
        u_face = self.u + u_dy * fn[1] 
        v_face = self.v + v_dy * fn[1] 
        p_face = self.p + p_dy * fn[1] 

        f.get_primitive_values(rho_face, u_face, v_face, p_face)


    def get_flux_and_apply(self, dt):


        for i, f in enumerate(self.faces):
            if i in (1,3):
                sign = -1
            else:
                sign = 1

            self.m -=  sign* f.surface * f.flux_Mass *dt
            self.mu -=  sign* f.surface * f.flux_Momx * dt
            self.mv -= sign*  f.surface * f.flux_Momy * dt
            self.e -=  sign* f.surface * f.flux_Energy * dt

    
        # update primitive values
        self.calc_primitives()

    def __str__(self):

        list_neighbors = ''
        for n in self.neighbors:
            list_neighbors = list_neighbors + f';{n.number}\n'
        return f"#Cell {self.number}, #center {self.center.__str__()}\n"+\
               f"with the #Neighbors {list_neighbors}\n"+\
               f"#Faces: {[f for f in self.faces]}"+\
               f'\n________________________________________________________________________'

    def _str__(self):
        non_ortho = np.array(self.non_orto_angles)*180/np.pi
        return f'#Cell {self.number}, \n'+\
                f'primitives: rho {self.rho}, u {self.u}, v {self.v}, p {self.p} \n'+\
                f'gradients: {self.gradients} \n'+\
                f'conserved: m, mu, mv, e {self.m, self.mu, self.mv, self.e}\n'+\
                f'non orthogonality to neighbors: {non_ortho}\n'+\
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

