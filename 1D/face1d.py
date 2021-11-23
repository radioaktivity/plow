import numpy as np
from global_proporties import *
from point import *

class Face1D:

    def __init__(self, number, center:Point) -> None:
        self.number = number
        self.center = center
        self.wormhole_face = None
        self.surface = None

        self.rho_L, self.u_L, self.p_L = None, None, None
        self.rho_R, self.u_R, self.p_R = None, None, None

        self.isL = False
        self.isR = False
        self.flux_Mass = None
        self.flux_Momx = None
        self.flux_Energy = None

    def getPrimitives(self, rho, u, p, side=None):
        # side: pov face so, when cell on the right side='R'

        if side == 'L':
            self.rho_L, self.u_L, self.p_L = rho, u, p
            self.isL = True
        else:
            self.rho_R, self.u_R, self.p_R = rho, u, p
            self.isR = True

    def calcFlux(self):
        w = self.wormhole_face
        if self.isL == False:
            self.rho_L, self.u_L, self.p_L = w.rho_L, w.u_L, w.p_L
        if self.isR == False:
            self.rho_R, self.u_R, self.p_R = w.rho_R, w.u_R, w.p_R
        
        self.flux_Mass, self.flux_Momx, self.flux_Energy = \
             self.getFlux(self.rho_L, self.rho_R, self.u_L, self.u_R, self.p_L, self.p_R)

        self.isL = False
        self.isR = False


    def getFlux(self, rho_L, rho_R, u_L, u_R, p_L, p_R):
        
        # left and right energies
        en_L = p_L/(atm.gamma-1)+0.5*rho_L * (u_L**2)
        en_R = p_R/(atm.gamma-1)+0.5*rho_R * (u_R**2)

        # compute star (averaged) states
        rho_star  = 0.5*(rho_L + rho_R)
        momx_star = 0.5*(rho_L * u_L + rho_R * u_R)
        en_star   = 0.5*(en_L + en_R)

        p_star = (atm.gamma-1)*(en_star-0.5*(momx_star**2)/rho_star)

        # compute fluxes (local Lax-Friedrichs/Rusanov)
        flux_Mass   = momx_star
        flux_Momx   = momx_star**2/rho_star + p_star
        flux_Energy = (en_star+p_star) * momx_star/rho_star

        # find wavespeeds
        C_L = np.sqrt(atm.gamma*p_L/rho_L) + np.abs(u_L)
        C_R = np.sqrt(atm.gamma*p_R/rho_R) + np.abs(u_R)
        C = np.maximum( C_L, C_R )

        # add stabilizing diffusive term
        flux_Mass   -= C * 0.5 * (rho_L - rho_R)
        flux_Momx   -= C * 0.5 * (rho_L * u_L - rho_R * u_R)
        flux_Energy -= C * 0.5 * ( en_L - en_R )

        return flux_Mass, flux_Momx, flux_Energy

    def __str__(self):
        return f'Face number {self.number} at %s'%(self.center)

if __name__ == '__main__':
    f1 = Face1D(1, Point(np.array([0,1])))
    print(f1)