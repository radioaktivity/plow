

def getConserved( rho, u, p, gamma, vol):

	Mass   = rho * vol
	Momx   = rho * u * vol
	Energy = (p/(gamma-1) + 0.5*rho*u**2)*vol

	return Mass, Momx, Energy


def getPrimitive( Mass, Momx, Energy, gamma, vol ):

	rho = Mass / vol
	u  = Momx / rho / vol
	p   = (Energy/vol - 0.5*rho * u**2) * (gamma-1)
	return rho, u, p
