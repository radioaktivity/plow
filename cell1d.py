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
        
    def calc_gradient(self):
        pass

    def extrapol_in_time(self, dt):
        # extrapolate half-step in time
        self.rho = self.rho - 0.5*dt * ( self.u * self.rho_dx + self.rho * self.u_dx)
        self.u  = self.u  - 0.5*dt * ( self.u * self.u_dx + (1/self.rho) * self.p_dx )
        self.p   = self.p   - 0.5*dt * ( atm.gamma* self.p * (self.u_dx)  + self.u * self.p_dx)
		

    def __str__(self):
        return f'cell number {self.number} between %s and %s'%(self.boundary_points[0], self.boundary_points[1])