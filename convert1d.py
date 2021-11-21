from global_proporties import *

def getConserved( rho, u, p, vol):
	
	Mass   = rho * vol
	Momx   = rho * u * vol
	Energy = (p/(atm.gamma-1) + 0.5*rho*u**2)*vol

	return Mass, Momx, Energy


def getPrimitive( Mass, Momx, Energy, vol ):

	rho = Mass / vol
	u  = Momx / rho / vol
	p   = (Energy/vol - 0.5*rho * u**2) * (atm.gamma-1)
	return rho, u, p
