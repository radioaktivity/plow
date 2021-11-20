import numpy as np

from face import *
from point import *
from convert1d import *

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
        self.volume = 0 # or length in 1D        

        # primitives
        self.rho = 0
        self.u = 0
        self.p = 0

        # gradients
        self.rho_dx = 0
        self.u_dx = 0
        self.p_dx = 0

        # conservatives
        self.Mass = 0
        self.Momx = 0
        self.Energy = 0

    def calc_primitives(self):
        self.rho, self.u, self.p = \
            getPrimitive(self.Mass, self.Momx, self.Energy, self.volume)
    
    def calc_center(self):
        self.center.Y = 0
        self.center.X = (self.boundary_points[0]+self.boundary_points[1])/2

    def calc_all(self):
        self.calc_center()

    def calc_gradient(self):

        d = 2* self.faces[2].surface 
        n_L = self.neighbors[0]
        n_R = self.neighbors[1]

        rho_dx = (n_R.rho-n_L.rho)/d 
        u_dx =  (n_R.u-n_L.u)/d
        v_dx =  (n_R.v-n_L.v)/d
        p_dx =  (n_R.p-n_L.p)/d


    def extrapol_in_time(self, dt):
        # extrapolate half-step in time
        self.rho = self.rho - 0.5*dt * ( self.u * self.rho_dx + self.rho * self.u_dx)
        self.u  = self.u  - 0.5*dt * ( self.u * self.u_dx + (1/self.rho) * self.p_dx )
        self.p   = self.p   - 0.5*dt * ( atm.gamma* self.p * (self.u_dx)  + self.u * self.p_dx)
		

    def __str__(self):
        return f'cell number {self.number} between %s and %s and neighbors {[n.number for n in self.neighbors]}'\
            %(self.boundary_points[0], self.boundary_points[1])