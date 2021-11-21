import numpy as np

from face import *
from point import *
from convert1d import *
import scipy.interpolate

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
        self.volume = None # or length in 1D        
        self.distance = None

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
    
    def calc_conserved(self):
        self.Mass, self. Momx, self.Energy = \
            getConserved(self.rho, self.u, self.p, self.volume)

    def calc_center(self):
        self.center.Y = 0
        self.center.X = (self.boundary_points[0].X+self.boundary_points[1].X)/2

    def calc_volume(self):
        self.volume = self.boundary_points[0].distance(self.boundary_points[1])


    def calc_all(self):
        self.calc_center()
        self.calc_volume()
        self.distance = self.volume # 1d case

    def calc_gradients(self, type='central'):
        
        if type == 'central':
            d = 2* self.distance
            n_L = self.neighbors[0]
            n_R = self.neighbors[1]

            self.rho_dx = (n_R.rho-n_L.rho)/d 
            self.u_dx =  (n_R.u-n_L.u)/d
            self.p_dx =  (n_R.p-n_L.p)/d
        elif type == 'forward':
            d = self.distance
            n_L = self.neighbors[0]
            n_R = self.neighbors[1]

            if self.u >=0:
                self.rho_dx = (self.rho-n_L.rho)/d 
                self.u_dx =  (self.u-n_L.u)/d
                self.p_dx =  (self.p-n_L.p)/d
            else:
                self.rho_dx = (self.rho-n_R.rho)/d 
                self.u_dx =  (self.u-n_R.u)/d
                self.p_dx =  (self.p-n_R.p)/d
        elif type == 'spline':
            n_L = self.neighbors[0]
            n_L_L = n_L.neighbors[0]
            n_R = self.neighbors[1]
            n_R_R = n_R.neighbors[1]

            #x = [n_L.center.X, n_L_L.center.X, self.center.X, n_R.center.X, n_R_R.center.X]
            x = [-self.volume*2, -self.volume,0,self.volume,2*self.volume]

            self.rho_dx = self.getGradientSpline(
                [n_L_L.rho, n_L.rho, self.rho, n_R.rho, n_R_R.rho],x)
            self.u_dx = self.getGradientSpline(
                [n_L_L.u, n_L.u, self.u, n_R.u, n_R_R.u],x)
            self.p_dx = self.getGradientSpline(
                [n_L_L.p, n_L.p, self.p, n_R.p, n_R_R.p],x)

    def getGradientSpline(self, f, x):         

        spl = scipy.interpolate.splrep(x,f,k=3) # no smoothing, 3rd order spline
        ddy = scipy.interpolate.splev(x,spl,der=1) # use those knots to get second derivative

        return ddy[2]


    def extrapol_in_time(self, dt):
        # extrapolate half-step in time
        self.rho = self.rho - 0.5*dt * ( self.u * self.rho_dx \
            + self.rho * self.u_dx)
        self.u  = self.u  - 0.5*dt * ( self.u * self.u_dx \
            + (1/self.rho) * self.p_dx )
        self.p   = self.p   - 0.5*dt * ( atm.gamma* self.p \
            * (self.u_dx)  + self.u * self.p_dx)

    def extrapol2faces(self):
        # Left Face
        f = self.faces[0]
        fn = self.center.getVecBetween(f.center)
        rho_face = self.rho + self.rho_dx * fn[0]
        u_face = self.u + self.u_dx * fn[0]
        p_face = self.p + self.p_dx * fn[0]
        f.getPrimitives(rho_face, u_face, p_face, side='R')

        # Right Face
        f = self.faces[1]
        fn = self.center.getVecBetween(f.center)
        rho_face = self.rho + self.rho_dx * fn[0]
        u_face = self.u + self.u_dx * fn[0]
        p_face = self.p + self.p_dx * fn[0]
        f.getPrimitives(rho_face, u_face, p_face, side='L')

    def applyFlux(self, dt):
        
        f_L = self.faces[0]
        f_R = self.faces[1]
        
        self.Mass -= dt * (f_L.flux_Mass - f_R.flux_Mass) * self.volume 
        self.Momx -= dt * (f_L.flux_Momx -f_R.flux_Momx) * self.volume 
        self.Energy -= dt * (f_L.flux_Energy - f_R.flux_Energy) * self.volume 

        self.calc_primitives()

    def print_values(self):
        string = f'cell {self.number} has rho: {self.rho} u: {self.u} and p: {self.p}'
        string += '\n____________________________________________'
        print(string)

    def __str__(self):
        return f'cell number {self.number} between %s and %s and neighbors {[n.number for n in self.neighbors]}'\
            %(self.boundary_points[0], self.boundary_points[1])