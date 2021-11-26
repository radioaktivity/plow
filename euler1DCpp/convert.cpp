
#include "convert.h"

float * getConserved(float rho,float u,float p,float vol){

  float gamma = 5/3;
	float Mass   = rho * vol;
	float Momx   = rho * u * vol;
	float Energy = (p/(gamma-1) + 0.5*rho*u*u)*vol;

  static float resultarray[3] = {Mass, Momx, Energy};

	return resultarray;
}

float * getPrimitive( float Mass, float Momx, float Energy, float vol ){
  float gamma = 5/3;
	float rho = Mass / vol;
	float u  = Momx / rho / vol;
	float p   = (Energy/vol - 0.5*rho * u*u) * (gamma-1);

  static float resultarray[3] = {rho, u, p};

	return resultarray;
}