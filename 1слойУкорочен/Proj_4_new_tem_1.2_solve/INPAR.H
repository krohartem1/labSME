// inpar.h
#ifndef _INPAR_LIBSAM

#define _INPAR_LIBSAM

// Definition of auxiliary structures
#include <stdio.h>

// Definition of auxiliary structures
class MaterialConstants{
 public:
  Real SY2; //yield limit for steel 
  Real SY1; // initial phase yield limit for TiNi at martensitic state
  Real SY1A; // initial phase yield limit for TiNi at austenitis state
  Real E2; // Young modulus of steel
  Real E1A; // Young_A
  Real E1M_L;//Young_M at loading
  Real E1M_U; // Young_M at unloading
  Real H2; // strengthening coefficient for steel
  Real H1; // strengthening coefficient for TiNi at martensitic state
  Real H1A; // strengthening coefficient for TiNi at austenitic state
  Real h1; // height of TiNi element
  Real h2; // height of TiNi element
  Real h; // total height h1+h2
  Real b1; // width of steel element
  Real b2; // width of TiNi element
  Real Alf1A; // coefficient of thermal expansion for TiNi element in mart. state
  Real Alf1M; // coefficient of thermal expansion for TiNi element in aust. state
  Real Alf2; // coefficient of thermal expansion for steel element

  Real k_b; // ratio b1/b2
  Real k_h; // ratio h1/(h1+h2)
  Real Kr; // recovery ratio
  Real EpsT; //preliminary strain
  Real r_l; // ratio gap to length
  Real l_n; // length of plate
  Real rad; // iniyial radius of curvature

  Real Mf; //
  Real Ms; //
  Real As; //
  Real Af; //

  Real Tmin; //minimum temperature of the experiment
  Real Tmax; //maximum temperature of the experiment
 
  Real Lam1; // coefficient for phase strain due to TWSM
  Real ETP; // coefficient for phase strain due to transformation plasticity

  Real M_t; // total moment per unit of depth, N 
//Default constructor
//  MaterialConstants(void){};
//Initialiser reading constants from file
  void Init(void);
//Default destructor
 // ~MaterialConstants(void){};
}; 
//end of class MaterialConstants
// ***************************************


#endif
