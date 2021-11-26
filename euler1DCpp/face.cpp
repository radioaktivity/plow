#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm> 
#include "point.h"
#include "face.h"


using namespace std;
Face::Face(){
  // Do nothing
}

Face::Face(int NUMBER, Point CENTER){
    number = NUMBER;
    center = CENTER;
}

void Face::getPrimitives(float rho, float u, float p, int side){
  isL = false;
  isR = false;
  if (side == 0){
      rho_L = rho;
      u_L = u;
      p_L = p;
      isL = true;
  }
  else{
      rho_R = rho;
      u_R = u;
      p_R = p;
      isR = true;
  }
}
void Face::calcFlux(){
  if (!isL){
    rho_L = wormhole_face->rho_L;
    u_L = wormhole_face->u_L;
    p_L = wormhole_face->p_L;
  }
  if (!isR){
    rho_R = wormhole_face->rho_R;
    u_R = wormhole_face->u_R;
    p_R = wormhole_face->p_R;
  }
  isL = false;
  isR = false;
  Face::getFlux();
}

void Face::getFlux(){
  float gamma = 5./3.;
  // left and right energies
  float en_L = p_L/( gamma-1)+0.5*rho_L * (u_L*u_L);
  float en_R = p_R/( gamma-1)+0.5*rho_R * (u_R*u_R);

  // compute star (averaged) states
  float rho_star  = 0.5*(rho_L + rho_R);
  float momx_star = 0.5*(rho_L * u_L + rho_R * u_R);
  float en_star   = 0.5*(en_L + en_R);

  float p_star = ( gamma-1)*(en_star-0.5*(momx_star*momx_star)/rho_star);

  // compute fluxes (local Lax-Friedrichs/Rusanov)
  flux_Mass   = momx_star;
  flux_Momx   = momx_star*momx_star/rho_star + p_star;
  flux_Energy = (en_star+p_star) * momx_star/rho_star;

  // find wavespeeds
  float C_L =sqrt( gamma*p_L/rho_L) + std::abs(u_L);
  float C_R = sqrt( gamma*p_R/rho_R) + std::abs(u_R);
  float C = std::max( C_L, C_R );

  // add stabilizing diffusive term
  flux_Mass  -= C * 0.5 * (rho_L - rho_R);
  flux_Momx   -= C * 0.5 * (rho_L * u_L - rho_R * u_R);
  flux_Energy -= C * 0.5 * ( en_L - en_R );
}
