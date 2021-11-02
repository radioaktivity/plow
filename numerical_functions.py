import numpy as np

def getGradient(f, dx):
    """
    Calculate the gradients of a field
    f        is a matrix of the field
    dx       is the cell size
    f_dx     is a matrix of derivative of f in the x-direction
    f_dy     is a matrix of derivative of f in the y-direction
    """
    # directions for np.roll() 
    R = -1   # right
    L = 1    # left

    f_dx = ( np.roll(f,R,axis=0) - np.roll(f,L,axis=0) ) / (2*dx)
    f_dy = ( np.roll(f,R,axis=1) - np.roll(f,L,axis=1) ) / (2*dx)

    return f_dx, f_dy

def extrapolateInSpaceToFace(f, f_dx, f_dy, dx):
    """
    Calculate the gradients of a field
    f        is a matrix of the field
    f_dx     is a matrix of the field x-derivatives
    f_dy     is a matrix of the field y-derivatives
    dx       is the cell size
    f_XL     is a matrix of spatial-extrapolated values on `left' face along x-axis 
    f_XR     is a matrix of spatial-extrapolated values on `right' face along x-axis 
    f_YR     is a matrix of spatial-extrapolated values on `left' face along y-axis 
    f_YR     is a matrix of spatial-extrapolated values on `right' face along y-axis 
    """
    # directions for np.roll() 
    R = -1   # right
    L = 1    # left

    f_XL = f - f_dx * dx/2
    f_XL = np.roll(f_XL,R,axis=0)
    f_XR = f + f_dx * dx/2

    f_YL = f - f_dy * dx/2
    f_YL = np.roll(f_YL,R,axis=1)
    f_YR = f + f_dy * dx/2

    return f_XL, f_XR, f_YL, f_YR

def getFlux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R, gamma):
    """
    Calculate fluxed between 2 states with local Lax-Friedrichs/Rusanov rule 
    rho_L        is a matrix of left-state  density
    rho_R        is a matrix of right-state density
    vx_L         is a matrix of left-state  x-velocity
    vx_R         is a matrix of right-state x-velocity
    vy_L         is a matrix of left-state  y-velocity
    vy_R         is a matrix of right-state y-velocity
    P_L          is a matrix of left-state  pressure
    P_R          is a matrix of right-state pressure
    gamma        is the ideal gas gamma
    flux_Mass    is the matrix of mass fluxes
    flux_Momx    is the matrix of x-momentum fluxes
    flux_Momy    is the matrix of y-momentum fluxes
    flux_Energy  is the matrix of energy fluxes
    """

    # left and right energies
    en_L = P_L/(gamma-1)+0.5*rho_L * (vx_L**2+vy_L**2)
    en_R = P_R/(gamma-1)+0.5*rho_R * (vx_R**2+vy_R**2)

    # compute star (averaged) states
    rho_star  = 0.5*(rho_L + rho_R)
    momx_star = 0.5*(rho_L * vx_L + rho_R * vx_R)
    momy_star = 0.5*(rho_L * vy_L + rho_R * vy_R)
    en_star   = 0.5*(en_L + en_R)

    P_star = (gamma-1)*(en_star-0.5*(momx_star**2+momy_star**2)/rho_star)

    # compute fluxes (local Lax-Friedrichs/Rusanov)
    flux_Mass   = momx_star
    flux_Momx   = momx_star**2/rho_star + P_star
    flux_Momy   = momx_star * momy_star/rho_star
    flux_Energy = (en_star+P_star) * momx_star/rho_star

    # find wavespeeds
    C_L = np.sqrt(gamma*P_L/rho_L) + np.abs(vx_L)
    C_R = np.sqrt(gamma*P_R/rho_R) + np.abs(vx_R)
    C = np.maximum( C_L, C_R )

    # add stabilizing diffusive term
    flux_Mass   -= C * 0.5 * (rho_L - rho_R)
    flux_Momx   -= C * 0.5 * (rho_L * vx_L - rho_R * vx_R)
    flux_Momy   -= C * 0.5 * (rho_L * vy_L - rho_R * vy_R)
    flux_Energy -= C * 0.5 * ( en_L - en_R )

    return flux_Mass, flux_Momx, flux_Momy, flux_Energy

