from point import *

class Face:
    def __init__(self, p1:Point, p2:Point):
        # geometry
        self.center = 0
        self.boundary_points = []
        self.sruface = 0

        # primitive values
        self.rho = 0
        self.p = 0
        self.u = 0
        self.v = 0
