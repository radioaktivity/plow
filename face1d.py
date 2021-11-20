import numpy as np

from point import *

class Face1D:

    def __init__(self, number, center:Point) -> None:
        self.number = number
        self.center = center
        self.wormhole_face = None

    def getFlux(rho_L, rho_R, u_L, u_R, p_L, p_R, gamma):
        # left and right energies
        en_L = p_L/(gamma-1)+0.5*rho_L * (u_L**2)
        en_R = p_R/(gamma-1)+0.5*rho_R * (u_R**2)

        # compute star (averaged) states
        rho_star  = 0.5*(rho_L + rho_R)
        momx_star = 0.5*(rho_L * u_L + rho_R * u_R)
        momy_star = 0.5*(rho_L * vy_L + rho_R * vy_R)
        en_star   = 0.5*(en_L + en_R)

        p_star = (gamma-1)*(en_star-0.5*(momx_star**2+momy_star**2)/rho_star)

        # compute fluxes (local Lax-Friedrichs/Rusanov)
        flux_Mass   = momx_star
        flux_Momx   = momx_star**2/rho_star + p_star
        flux_Energy = (en_star+p_star) * momx_star/rho_star

        # find wavespeeds
        C_L = np.sqrt(gamma*p_L/rho_L) + np.abs(u_L)
        C_R = np.sqrt(gamma*p_R/rho_R) + np.abs(u_R)
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