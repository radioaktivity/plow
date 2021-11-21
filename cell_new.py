from matplotlib.pyplot import delaxes
import numpy as np
import random

from numpy import linalg

from point import *
from face_new import *
import convert
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
        
        self.rho_dx, self.rho_dy = 0, 0
        self.u_dx, self.u_dy = 0, 0
        self.v_x, self.v_dy = 0, 0
        self.p_dx, self.p_dy = 0, 0

        # conservatives
        self.Mass = 0
        self.Momx = 0
        self.Momy = 0
        self.Energy = 0

    ############################
    ## Functions for meshing ###
    ############################

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
        if(len(self.boundary_points)==3):
            raise Exception('Used triangle')
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

            self.longest_side = min([f.surface for f in self.faces])
            self.volume = surface  

    def calc_vecs2faces(self):
        # Calculates distances to all face centers of the cell
        
        for f in self.faces:
            self.vecs2faces.append(self.center.getVecBetween(f.center))

    ##############################
    #### Functions for Solver ####
    ##############################

    def calc_all(self):
        self.calc_volume()
        self.calc_vecs2faces()

    def calc_primitives(self):
        self.rho, self.u, self.v, self.p = \
            convert.getPrimitive(self.Mass, self.Momx, self.Momy, self.Energy, self.volume)

    def calc_conserved(self):
        self.Mass, self.Momx, self.Momy, self.Energy = \
            convert.getConserved(self.rho, self.u, self.v, self.p, self.volume)

    def calc_gradients_central(self):

        n_L = self.neighbors[0] 
        n_R = self.neighbors[1]

        # top face length 
        #d = -1 * self.center.distance(n_L.center) + self.center.distance(n_R.center)
        d = 2*self.longest_side

        self.rho_dx = (n_R.rho-n_L.rho)/d 
        self.u_dx =  (n_R.u-n_L.u)/d
        self.v_dx =  (n_R.v-n_L.v)/d
        self.p_dx =  (n_R.p-n_L.p)/d

        n_D = self.neighbors[3]
        n_U = self.neighbors[2]

        # left face length
        # d = -1*  self.center.distance(n_U.center) + self.center.distance(n_D.center)
        d = 2*self.longest_side

        self.rho_dy = (n_U.rho-n_D.rho)/d
        self.u_dy =  (n_U.u-n_D.u)/d
        self.v_dy =  (n_U.v-n_D.v)/d
        self.p_dy = (n_U.p-n_D.p)/d
 
    def extrapol_in_time(self, dt):

        # extrapolate half-step in time
        self.rho = self.rho - 0.5*dt * \
            ( self.u * self.rho_dx + self.rho * self.u_dx \
                + self.v * self.rho_dy + self.rho * self.v_dy)
        self.u = self.u - 0.5*dt * \
            ( self.u * self.u_dx + self.v * self.u_dy \
                + (1/self.rho) * self.p_dx )
        self.v = self.v - 0.5*dt * \
            ( self.u * self.v_dx + self.v * self.v_dy \
                + (1/self.rho) * self.p_dy )
        self.p = self.p - 0.5*dt * \
            ( atm.gamma*self.p * (self.u_dx + self.v_dy)\
                  + self.u * self.p_dx + self.v * self.p_dy )

    def extrapol2faces(self):

        #left right top bottom
        # x x y y 
        # LEFT
        f = self.faces[0]

        fn = self.vecs2faces[0]

        rho_face = self.rho + self.rho_dx * fn[0]
        u_face = self.u + self.u_dx * fn[0]
        v_face = self.v + self.v_dx * fn[0]
        p_face = self.p + self.p_dx * fn[0]
            
        f.get_primitive_values(rho_face, u_face, v_face, p_face, cellsent='R')

        # RIGHT
        f = self.faces[1]

        fn = self.vecs2faces[1]

        rho_face = self.rho + self.rho_dx * fn[0]
        u_face = self.u + self.u_dx * fn[0]
        v_face = self.v + self.v_dx * fn[0]
        p_face = self.p + self.p_dx * fn[0]

        f.get_primitive_values(rho_face, u_face, v_face, p_face, cellsent='L')

        # TOP
        f = self.faces[2]

        fn = self.vecs2faces[2]

        rho_face = self.rho + self.rho_dy * fn[1]
        u_face = self.u + self.u_dy * fn[1]
        v_face = self.v + self.v_dy * fn[1]
        p_face = self.p + self.p_dy * fn[1]

        f.get_primitive_values(rho_face, u_face, v_face, p_face, cellsent='BOTTOM')


        # BOTTOM
        f = self.faces[3]

        fn = self.vecs2faces[3]

        rho_face = self.rho + self.rho_dy * fn[1]
        u_face = self.u + self.u_dy * fn[1]
        v_face = self.v + self.v_dy * fn[1]
        p_face = self.p + self.p_dy * fn[1]

        f.get_primitive_values(rho_face, u_face, v_face, p_face, cellsent='TOP')


    def get_flux_and_apply(self, dt):
        # L R U D
        f_L = self.faces[0]
        f_R = self.faces[1]
        self.Mass -= dt * (f_L.flux_Mass - f_R.flux_Mass) * f_L.surface 
        self.Momx -= dt * (f_L.flux_Momx -f_R.flux_Momx) * f_L.surface 
        self.Momy -= dt * (f_L.flux_Momy -f_R.flux_Momy) * f_L.surface  
        self.Energy -= dt * (f_L.flux_Energy - f_R.flux_Energy) * f_L.surface  
     
        f_L = self.faces[2]
        f_R = self.faces[3]
        self.Mass -= dt * (f_L.flux_Mass - f_R.flux_Mass) * f_L.surface  
        self.Momx -= dt * (f_L.flux_Momx -f_R.flux_Momx) * f_L.surface  
        self.Momx -= dt * (f_L.flux_Momy -f_R.flux_Momy) * f_L.surface  
        self.Energy -= dt * (f_L.flux_Energy - f_R.flux_Energy) * f_L.surface  
     

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
                f'conserved: m, mu, mv, e {self.Mass, self.Momx, self.Momy, self.Energy}\n'+\
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

