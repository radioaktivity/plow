import matplotlib.pyplot as plt
from convert import getConserved

from point import *
from global_proporties import *
from vector_alg import *

class Face:
    def __init__(self, p1:Point, p2:Point):
        if p1 is p2:
            raise Exception("Face generation with 2 equal points")

        # geometry
        self.center = Point()
        self.boundary_points = [p1, p2]
        self.calc_center()
        self.surface = p1.distance(p2)
        self.n = [0, 0]
        self.theta = 0 # angle to x axis
        self.calc_surfacenormal() 
        
        # cells conected
        self.cells_connected = []

        # primitive values
        [self.rho_L_X, self.rho_R_X, self.u_L_X, self.u_R_X, self.v_R_X, self.v_L_X, self.p_L_X, self.p_R_X]=\
            [0.,0.,0.,0.,0.,0.,0.,0.]
        [self.rho_L_Y, self.rho_R_Y, self.u_L_Y, self.u_R_Y, self.v_R_Y, self.v_L_Y, self.p_L_Y, self.p_R_Y]=\
            [0.,0.,0.,0.,0.,0.,0.,0.]
        
        self.flux_Mass_X  = 0
        self.flux_Momx_X   = 0
        self.flux_Momy_X   = 0
        self.flux_Energy_X = 0
        self.flux_Mass_Y  = 0
        self.flux_Momx_Y   = 0
        self.flux_Momy_Y   = 0
        self.flux_Energy_Y = 0

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
        self.cells_connected.append(cell)
    
    def get_neighbor(self, cell):
        '''
        returns the neighbor cell

        self.cells_connected contains two cells
        the cell that sends the demand 
        and the cell which is on the other side 
        of the face

        the cell which sends the demand wants the cell object on the other side returned
        '''
        if self.cells_connected[0] == cell:
            if len(self.cells_connected) == 2:
                return self.cells_connected[1]
            else:
                return False # no neighbor
        else:
            return self.cells_connected[0]

    def calc_surfacenormal(self):
        # Calculate a vector normal to the face 
        tangent = self.boundary_points[1].getVec()-self.boundary_points[0].getVec()

        # normal of [a,b] -> 1/norm([a,b]) * [-b, a]
        normal = [-1/self.surface * tangent[1], 
                    1/self.surface *  tangent[0]]
        self.n = normal
        self.theta = angle_between( np.array([0,1]) , self.n )

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
    
    def get_primitive_values(self, rho_face_X, u_face_X, v_face_X, p_face_X,
                                    rho_face_Y, u_face_Y, v_face_Y, p_face_Y):
        # takes primitive values from one cell 
        # assigns primitive values to L 
        # if L is already full assigns to R

        if not(self.isL):
            self.rho_L_X, self.rho_L_Y = rho_face_X, rho_face_Y
            self.u_L_X, self.u_L_Y = u_face_X, u_face_Y
            self.v_L_X, self.v_L_Y = v_face_X, v_face_Y
            self.p_L_X, self.p_L_Y = p_face_X, p_face_Y

            self.isL = True
        else:
            self.rho_R_X, self.rho_R_Y = rho_face_X, rho_face_Y
            self.u_R_X, self.u_R_Y = u_face_X, u_face_Y
            self.v_R_X, self.v_R_Y = v_face_X, v_face_Y
            self.p_R_X, self.p_R_Y = p_face_X, p_face_Y

            self.isR = True

    def getFlux(self, gamma = 5/3):
        if self.isR: # is inner face
            self.flux_Mass_X, \
                self.flux_Momx_X, \
                self.flux_Momy_X, \
                self.flux_Energy_X = \
                    self.calcFlux(self.rho_L_X, self.rho_R_X, 
                                self.u_L_X, self.u_R_X, 
                                self.v_L_X, self.v_R_X, 
                                self.p_L_X, self.p_R_X)

            # notice here momy and momx changed place as well as v_L_Y,v_R_Y and u_L_Y,u_R_Y 
            self.flux_Mass_Y, \
                self.flux_Momy_Y, \
                self.flux_Momx_Y, \
                self.flux_Energy_Y = \
                    self.calcFlux(self.rho_L_Y, self.rho_R_Y,
                                  self.v_L_Y, self.v_R_Y,
                                  self.u_L_Y, self.u_R_Y, 
                                  self.p_L_Y, self.p_R_Y)

            self.isR = False
            self.isL = False # Reset the state so next iteration overwrites 
        else: # is boundary face
            # Don't change anything
            # Their standard value is 0
            self.isL = False # Reset the state so next iteration overwrites 
            
    def calcFlux(self, rho_L, rho_R, u_L, u_R, v_R, v_L, p_L, p_R ):
        gamma = atm.gamma
        if (p_L*p_R)<0:
            raise Exception("negative pressure")
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

    def display_normals(self):
        display_vector(self.center.getVec(), self.n, color='k')

    
    def calcFlux2(self, gamma):
        [rho_L, rho_R, u_L, u_R, v_R, v_L, p_L, p_R]= \
            [self.rho_L, self.rho_R, self.u_L, self.u_R, self.v_R, self.v_L, self.p_L, self.p_R]

        m_L, mu_L, mv_L, e_L = getConserved(rho_L,u_L,v_L,p_L)
        m_R, mu_R, mv_R, e_R = getConserved(rho_R,u_R,v_R,p_R)

        # find wavespeeds
        C_L = np.sqrt(gamma*p_L/rho_L) + np.abs(u_L)
        C_R = np.sqrt(gamma*p_R/rho_R) + np.abs(u_R)
        C = np.maximum( C_L, C_R )


        # F = 0.5 (F_L+F_R) - 0.5 * C * (u_R-u_L)
        flux_Mass   = 0.5 * (m_L+m_R) - C * 0.5 * (rho_L - rho_R)
        flux_Momx   = 0.5 *(mu_L+mu_R) - C * 0.5 * (rho_L * u_L - rho_R * u_R)
        flux_Momy   = 0.5 * (mv_L+mv_R)- C * 0.5 * (rho_L * v_L - rho_R * v_R)
        flux_Energy = 0.5 * (e_L+e_R) - C * 0.5 * ( e_L - e_R )

        return flux_Mass, flux_Momx, flux_Momy, flux_Energy

    


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