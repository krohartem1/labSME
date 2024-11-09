#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defglob.h"
#include "diaproc.h"
#include "inpar.h"
#define MALBUF 16
#define MAXBUF 128
// ------------------new-------------------------------------
const Real pi=3.14159265359;

// _____________MaterialConstants::Init________________________
void MaterialConstants::Init(void)
{
  FILE *fc;
  Real Dx;
  if( (fc=fopen("m_const.d","r"))==NULL){
    fprintf(stderr,"Cannot open m_const.d \n");
    pause();
	exit(1);
  }

  fkeyread_double(fc, "SY1", &Dx);    SY1 = Dx * 1E6; //from MPa into Pa
  fkeyread_double(fc, "SY2", &Dx);    SY2 = Dx * 1E6; //from MPa into Pa
  fkeyread_double(fc, "SY1A", &Dx);    SY1A = Dx * 1E6; //from MPa into Pa
  fkeyread_double(fc, "E2", &Dx); E2 = Dx * 1E9; //from GPa into Pa
  fkeyread_double(fc, "E1A", &Dx); E1A = Dx * 1E9; //from GPa into Pa
  fkeyread_double(fc, "E1M_L", &Dx); E1M_L = Dx * 1E9; //from GPa into Pa
  fkeyread_double(fc, "E1M_U", &Dx); E1M_U = Dx * 1E9; //from GPa into Pa
  fkeyread_double(fc, "H1", &Dx); H1 = Dx * 1E9; //from GPa into Pa 
  fkeyread_double(fc, "H1A", &Dx); H1A = Dx * 1E9; //from GPa into Pa 
  fkeyread_double(fc, "H2", &Dx); H2 = Dx * 1E9; //from GPa into Pa
//fkeyread_double(fc, "h", &Dx); h = Dx * 1E-3; //from mm into m
//fkeyread_double(fc, "b1", &Dx); b1 = Dx * 1E-3; //from mm into m
//fkeyread_double(fc, "b2", &Dx); b2 = Dx * 1E-3; //from mm into m
  fkeyread_double(fc, "h1", &Dx); h1 = Dx * 1E-3; //from mm into m
  fkeyread_double(fc, "h2", &Dx); h2 = Dx * 1E-3; //from mm into m
  fkeyread_double(fc, "k_h", &k_h); // ratio h1/(h1+h2)
  fkeyread_double(fc, "k_b", &k_b);
  fkeyread_double(fc, "Kr", &Kr); 
  fkeyread_double(fc, "EpsT", &Dx); EpsT = Dx * 1e-2; // from pst intu units
  fkeyread_double(fc, "r_l", &Dx); r_l = Dx * 1E-3;
  fkeyread_double(fc, "l_n", &Dx); l_n = Dx * 1E-3; //from mm into m
  fkeyread_double(fc, "rad", &Dx); rad = Dx * 1E-3; //from mm into m
  fkeyread_double(fc, "Mf", &Mf); 
  fkeyread_double(fc, "Ms", &Ms);
  fkeyread_double(fc, "As", &As);
  fkeyread_double(fc, "Af", &Af);
  fkeyread_double(fc, "Tmin", &Tmin); //minimum temperature of the experiment
  fkeyread_double(fc, "Tmax", &Tmax); //maximum temperature of the experiment
  fkeyread_double(fc, "Lam1", &Lam1); // coefficient for phase strain due to TWSM
  fkeyread_double(fc, "ETP", &Dx); ETP = Dx * 1E9;// coefficient for phase strain due to transformation plasticity from MPa into Pa
  fkeyread_double(fc, "M_t", &M_t);  //total moment per unit of depth, N 
  fkeyread_double(fc, "Alf1A", &Dx); Alf1A = Dx * 1e-6;//Heat expansion of TiNi in mart.state, 10^-6 K^-1
  fkeyread_double(fc, "Alf1M", &Dx); Alf1M = Dx * 1e-6;//Heat expansion of TiNi in aust.state, 10^-6 K^-1
  fkeyread_double(fc, "Alf2", &Dx); Alf2 = Dx * 1e-6;//Heat expansion of TiNi,  10^-6 K^-1

//    h1 = k_h*h;
//    h2 = h*(1.0-k_h);
	h = h1+h2;
	k_h = h1/h;


  printf("\n\nM_const INIT:");
  printf("\nSY1=%lg\tSY2=%lg\tSY1A=%lg"
		   ,SY1*1e-6,SY2*1e-6, SY1A*1e-6);
  printf("\nE1A=%9.5lg GPa\tE1M_L=%9.5lg GPa\tE1M_U=%9.5lg GPa\tE2=%9.5lg GPa" 
	       ,E1A*1e-9       ,E1M_L*1e-9       ,E1M_U*1e-9		,E2*1e-9);
  printf("\nH1=%9.5lg GPa\tH2=%9.5lg GPa\tH1A=%9.5lg GPa" 
	       ,H1*1e-9       ,H2*1e-9,      H1A*1e-9);
  printf("\nk_h=%lg" ,k_h);
// b1=b2*k_b;
  printf("\th1=%lg mm\th2=%lg mm\th=%lg mm" ,h1*1e3,h2*1e3,h*1e3);
//  printf("\nEpsT=%lg %", EpsT*1e2);
  printf("\nLam1=%lg\tETP=%9.5lg GPa"
           ,Lam1	,ETP*1e-9);
  printf("\nKr=%lg" ,Kr);
  printf("\n\nMs=%lg K\tMf=%lg K\tAs=%lg K\tAf=%lg K" ,Ms ,Mf ,As ,Af);


 	

   //pause();
//=================================================================
//=================================================================
  

} //end of Init for MaterialConstants




