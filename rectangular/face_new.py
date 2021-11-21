import matplotlib.pyplot as plt
from convert import getConserved

from point import *
from global_proporties import *
from vector_alg import *

class Face:
    def __init__(self, p1:Point, p2:Point, type='X'):
        if p1 is p2:
            raise Exception("Face generation with 2 equal points")
        
        # organisation
        self.number = 0

        # geometry
        self.center = Point()
        self.boundary_points = [p1, p2]
        self.calc_center()
        self.surface = abs(p1.distance(p2))
        self.n = [0, 0]
        self.theta = 0 # angle to x axis
        self.faceType = type
        self.wormhole_face = None

        self.calc_surfacenormal2() 

        # cells conected
        self.cells_connected = set()

        # primitive values
        [self.rho_L , self.rho_R , self.u_L , self.u_R , self.v_R , self.v_L , self.p_L , self.p_R ]=\
            [None]*8
        
        self.flux_Mass   = 0
        self.flux_Momx    = 0
        self.flux_Momy    = 0
        self.flux_Energy  = 0

        # logic switches
        # When 1. cell hands his values to this face:
        # isL is set to True
        # When 2. cell hands his values to this face
        # isR is also set to True
        # After the Flux calculation both variables are reset to zero
        self.isL = False
        self.isR = False

    def on_cell(self, cell):
        # adds a cell to the cells_connected array to later know
        # which cells this face is concading
        self.cells_connected.add(cell)

    def calc_surfacenormal(self):
        # Calculate a vector normal to the face 
        tangent = self.boundary_points[1].getVec()-self.boundary_points[0].getVec()

        # normal of [a,b] -> 1/norm([a,b]) * [-b, a]
        normal = [-1/self.surface * tangent[1], 
                    1/self.surface *  tangent[0]]
        self.n = normal
        self.theta = angle_between( np.array([0,1]) , self.n )
    
    def calc_surfacenormal2(self):
        if (self.faceType == 'X'):
            self.n = np.array([1,0])
        else:
            self.n = np.array([0,1])

    def calc_center(self):
        # Takes boundary poitns and sets center 
        boundary_points_ko =np.zeros([len(self.boundary_points),2])
        k = 0
        for p in self.boundary_points:
            boundary_points_ko[k] = p.getValue()
            k+=1
        self.center  = Point(boundary_points_ko.mean(axis=0))

    def is_equal_to(self, other_face):
        # Takes another Face() as input and checks its points
        # if the points are equal the face declares itself as equal
        i = 0 # i = 0 not identical; i = 1 one point identical; i = 2 both points identical 
        p1 = self.boundary_points[0]
        p2 = self.boundary_points[1]
        p1o = other_face.boundary_points[0]
        p2o = other_face.boundary_points[1]
        
        if p1 == p1o:
            i +=1
        if p1 == p2o:
            i +=1
        if p2 == p1o:
            i +=1
        if p2 == p2o:
            i +=1
        
        if i == 2:
            return True
        else:
            return False
    
    def get_primitive_values(self, rho, u, v, p, cellsent=None):
        # takes primitive values from one cell 
        # assigns primitive values to L 
        # if L is already full assigns to R

        if (cellsent == 'L') or (cellsent == 'BOTTOM'):
            self.rho_L = rho
            self.u_L = u
            self.v_L = v
            self.p_L = p

            self.isL = True

        else: # (cellsent == 'R') or (cellsent == 'TOP')
            self.rho_R = rho
            self.u_R = u
            self.v_R = v
            self.p_R = p

            self.isR = True

    def getFlux(self):


        if self.isR == False:
            self.rho_R = self.wormhole_face.rho_R
            self.u_R = self.wormhole_face.u_R
            self.v_R = self.wormhole_face.v_R
            self.p_R = self.wormhole_face.p_R
            self.isR = True
        elif self.isL == False:
            self.rho_L = self.wormhole_face.rho_L
            self.u_L = self.wormhole_face.u_L
            self.v_L = self.wormhole_face.v_L
            self.p_L = self.wormhole_face.p_L            
            self.isL = True

        if (self.isR and self.isL): # is inner face
            if (self.faceType == 'X'):
                self.flux_Mass, \
                    self.flux_Momx, \
                    self.flux_Momy, \
                    self.flux_Energy = \
                        self.calcFlux(self.rho_L, self.rho_R, 
                                    self.u_L, self.u_R, 
                                    self.v_L, self.v_R, 
                                    self.p_L, self.p_R)
            else: # self.faceType == 'Y'
                # notice here momy and momx changed place as well as v_L_Y,v_R_Y and u_L_Y,u_R_Y 
                self.flux_Mass, \
                    self.flux_Momy, \
                    self.flux_Momx, \
                    self.flux_Energy = \
                        self.calcFlux(self.rho_L, self.rho_R,
                                    self.v_L, self.v_R,
                                    self.u_L, self.u_R, 
                                    self.p_L, self.p_R)

            self.isR = False
            self.isL = False # Reset the state so next iteration overwrites 

    def calcFlux(self, rho_L, rho_R, u_L, u_R, v_R, v_L, p_L, p_R ):
        gamma = atm.gamma

        # left and right energies
        en_L = p_L/(gamma-1)+0.5*rho_L * (u_L**2+v_L**2)
        en_R = p_R/(gamma-1)+0.5*rho_R * (u_R**2+v_R**2)

        # compute star (averaged) states
        rho_star  = 0.5*(rho_L + rho_R)
        momx_star = 0.5*(rho_L * u_L + rho_R * u_R)
        momy_star = 0.5*(rho_L * v_L + rho_R * v_R)
        en_star   = 0.5*(en_L + en_R)
        
        P_star = (gamma-1)*(en_star-0.5*(momx_star**2+momy_star**2)/rho_star)
        
        # compute fluxes (local Lax-Friedrichs/Rusanov)
        flux_Mass   = momx_star
        flux_Momx   = momx_star**2/rho_star + P_star
        flux_Momy   = momx_star * momy_star/rho_star
        flux_Energy = (en_star+P_star) * momx_star/rho_star
        
        # find wavespeeds
        C_L = np.sqrt(gamma*p_L/rho_L) + np.abs(u_L)
        C_R = np.sqrt(gamma*p_R/rho_R) + np.abs(u_R)
        C = np.maximum( C_L, C_R )
        # plt.text(self.center.X, self.center.Y-atm.linespacing*1, f'WaveSpeed:{np.round(C, decimals=2)}')
        
        # add stabilizing diffusive term
        flux_Mass   -= C * 0.5 * (rho_L - rho_R)
        flux_Momx   -= C * 0.5 * (rho_L * u_L - rho_R * u_R)
        flux_Momy   -= C * 0.5 * (rho_L * v_L - rho_R * v_R)
        flux_Energy -= C * 0.5 * ( en_L - en_R )

        return flux_Mass, flux_Momx, flux_Momy, flux_Energy

    def give_FLuxes(self, side=None):
        if (side == 'L') or (side == 'TOP'):
            return  -self.flux_Mass, \
                        -self.flux_Momx, \
                        -self.flux_Momy, \
                        -self.flux_Energy
        else:
            return  self.flux_Mass, \
                        self.flux_Momx, \
                        self.flux_Momy, \
                        self.flux_Energy

    def display_normals(self):
        display_vector(self.center.getVec(), self.n, color='k')


    def plot_border(self):
        plt.plot([self.boundary_points[0].X, self.boundary_points[1].X], 
                [self.boundary_points[0].Y, self.boundary_points[1].Y], 'y')


    def __str__(self):
        return f'Face {str(hex(id(self)))}\n'+\
                f'Cells connected {self.cells_connected}\n'+\
                '____________________________________________________________'

    def _str__(self):
        return f"Face {str(hex(id(self)))}\n"+\
            '[rho_L, rho_R, u_L, u_R, v_R, v_L, p_L, p_R]= \n'+\
            f'{ [self.rho_L, self.rho_R, self.u_L, self.u_R, self.v_R, self.v_L, self.p_L, self.p_R]}\n'+\
                '___________________________________________'


if __name__ == "__main__":
    p1 = Point(np.array([1,1]))
    p2 = Point(np.array([3,2]))
    f1 = Face(p1,p2)
    f2 = Face(p2,p1)
    print(f1.is_equal_to(f2))
    print(f1.n)