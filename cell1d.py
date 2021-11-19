import numpy as np

from face import *
from point import *

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
        self.gradients = None

        # conservatives
        self.Mass = 0
        self.Momx = 0
        self.Energy = 0

    def calc_primitives(self):
        self.rho, self.u, self.p = \
            getPrimitive(self.Mass, self.Momx, self.Energy, self.volume)
        
    def __str__(self):
        return f'cell number {self.number} between %s and %s'%(self.boundary_points[0], self.boundary_points[1])