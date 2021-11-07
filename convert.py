def getConserved( rho, u, v, p, gamma, vol ):
    """
    Calculate the conserved variable from the primitive
    rho      is matrix of cell densities
    u       is matrix of cell x-velocity
    v       is matrix of cell y-velocity
    p        is matrix of cell pressures
    gamma    is ideal gas gamma
    vol      is cell volume
    m     is matrix of m in cells
    mu     is matrix of x-momentum in cells
    mv     is matrix of y-momentum in cells
    e   is matrix of e in cells
    """
    m   = rho * vol
    mu   = rho * u * vol
    mv   = rho * v * vol
    e = (p/(gamma-1) + 0.5*rho*(u**2+v**2))*vol

    return m, mu, mv, e

def getPrimitive( m, mu, mv, e, gamma, vol ):
    """
    Calculate the primitive variable from the conservative
    m     is matrix of m in cells
    mu     is matrix of x-momentum in cells
    mv     is matrix of y-momentum in cells
    e   is matrix of e in cells
    gamma    is ideal gas gamma
    vol      is cell volume
    rho      is matrix of cell densities
    u       is matrix of cell x-velocity
    v       is matrix of cell y-velocity
    p        is matrix of cell pressures
    """
    rho = m / vol
    u  = mu / rho / vol
    v  = mv / rho / vol
    p   = (e/vol - 0.5*rho * (u**2+v**2)) * (gamma-1)

    return rho, u, v, p