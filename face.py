from numpy import true_divide
from point import *
from global_proporties import *

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
        self.calc_surfacenormal()
        
        # cells conected
        self.cells_connected = []

        # primitive values
        [self.rho_L, self.rho_R, self.u_L, self.u_R, self.v_R, self.v_L, self.p_L, self.p_R]=\
            [0.,0.,0.,0.,0.,0.,0.,0.]
        
        # conserved values
        [self.m, self.mu, self.mv, self.e] = \
            [0.,0.,0.,0.]

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
    
    def get_vol_of_neighbor(self, cell):
        '''
        returns the volume of the neighbor cell

        self.cells_connected contains two cells
        the cell that sends the demand 
        and the cell which is on the other side 
        of the face
        '''
        if self.cells_connected[0] == cell:
            return self.cells_connected[1].volume
        else:
            return self.cells_connected[0].volume

    def calc_surfacenormal(self):
        # Calculate a vector normal to the face 
        tangent = self.boundary_points[1].getVec()-self.boundary_points[0].getVec()
        normal = [-1/self.surface * tangent[1], 
                    1/self.surface *  tangent[0]]
        self.n = normal

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
    
    def get_primitive_value(self, rho, u, v, p):
        # takes primitive values from one cell 
        # assigns primitive values to L 
        # if L is already full assigns to R

        if not(self.isL):
            self.rho_L = rho
            self.u_L = u 
            self.v_L = v
            self.p_L = p

            self.isL = True
        else:
            self.rho_R = rho
            self.u_R = u
            self.v_R = v
            self.p_R = p
            
            self.isR = True

    def getFlux(self, gamma = 5/3):
        if self.isR: # is inner face
            self.calcFlux(gamma)
            self.isR = False
            self.isL = False # Reset the state so next iteration overwrites 
        else: # is boundary face
            [self.m, self.mu, self.mv, self.e] = [0,0,0,0]
            self.isL = False # Reset the state so next iteration overwrites 
        
    
    def calcFlux(self, gamma):
        [rho_L, rho_R, u_L, u_R, v_R, v_L, p_L, p_R]= \
         [self.rho_L, self.rho_R, self.u_L, self.u_R, self.v_R, self.v_L, self.p_L, self.p_R]

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
        m   = momx_star
        mu   = momx_star**2/rho_star + P_star
        mv   = momx_star * momy_star/rho_star
        e = (en_star+P_star) * momx_star/rho_star


        # find wavespeeds
        C_L = np.sqrt(gamma*p_L/rho_L) + np.abs(u_L)
        C_R = np.sqrt(gamma*p_R/rho_R) + np.abs(u_R)
        C = np.maximum( C_L, C_R )

        # add stabilizing diffusive term
        m   -= C * 0.5 * (rho_L - rho_R)
        mu   -= C * 0.5 * (rho_L * u_L - rho_R * u_R)
        mv   -= C * 0.5 * (rho_L * v_L - rho_R * v_R)
        e -= C * 0.5 * ( en_L - en_R )
        
        [self.m, self.mu, self.mv, self.e] = \
        [m, mu, mv, e]

    def __str__(self):
        return f"Face between {self.boundary_points[0].__str__()} and {self.boundary_points[1].__str__()} center : {self.center}"


if __name__ == "__main__":
    p1 = Point(np.array([1,1]))
    p2 = Point(np.array([3,2]))
    f1 = Face(p1,p2)
    f2 = Face(p2,p1)
    print(f1.is_equal_to(f2))
    print(f1.n)