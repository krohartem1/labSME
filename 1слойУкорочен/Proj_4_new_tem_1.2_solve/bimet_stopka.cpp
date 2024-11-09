#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "stdafx.h"
//#include <winuser.h>
#include <windows.h>
//#include <iostream>
#include <dos.h>
#include <math.h>
//#include <conio.h>     // getch
#include "defglob.h"
#include "diaproc.h"
#include "inpar.h"

		int main ()
		{
			
			FILE *res; 
			FILE *res_sm;
			MaterialConstants MC;
			MC.Init(); //read constants from file m_const.d

			double Tem, hTem, Fi, hFi;
			double Eps1out, Eps1in, Eps2in, Eps2out, Eps1outF, Eps1inF, Eps1outMP, Eps1inMP; //Eps1E, Eps1F, hEps1F, Eps1P, hEps1P;
	//		double hEps1out, hEps1in, hEps2in, hEps2out;
			double Pt; //point on Fi-T curve
	//		double D1_old, D1_new, hD1; // compliance of TiNi

			double M, hM, Sig1, S1, Sig2, S2, hSig1, hS1, hSig2, hS2; 
			double Sig1out_new,	Sig1in_new,	Sig2in_new,	Sig2out_new;
			//double SigElastic; // stress in elastic approach SigElastic = 6*M/h**2
			double EpsElastic; // strain in elastic approach EpsElastic = h/2/R; 

			double Hev1out_old, Hev1in_old, Hev2out_old, Hev2in_old;
			double Hev1out_new, Hev1in_new, Hev2out_new, Hev2in_new;
			double Hev1out, Hev1in, Hev2out, Hev2in;

			double Podat1L = 1.0/MC.H1-1.0/MC.E1M_L; // reduced elastic-unelastic compliance of TiNi on loading
			double Podat1U = 1.0/MC.H1-1.0/MC.E1M_U; // reduced elastic-unelastic compliance of TiNi on unloading
			double Podat1A = 1.0/MC.H1A - 1.0/MC.E1A; // reduced elastic-unelastic compliance of TiNi at austenitic state
			double Podat2 = 1.0/MC.H2-1.0/MC.E2; // reduced elastic-unelastic compliance of steel
			double D2 = 1.0/MC.E2;
			double D1M_L = 1.0/MC.E1M_L;
			double D1M_U = 1.0/MC.E1M_U;
			double D1A = 1.0/MC.E1A;
			double P1 = 1.0/MC.H1;
			double P2 = 1.0/MC.H2;


			double TWSM = 0.0;
			double Eps1Fin_MF = 0.0;
			double Eps1Fin_MS = 0.0;
			double Eps1Fout_MF = 0.0;
			double Eps1Fout_MS = 0.0;

			double Eps0outAf = 0;
			double Eps0inAf = 0;

			double hEps1out = 0.0;
			double hEps2out = 0.0;
			double hEps1in = 0.0;
			double hEps2in = 0.0;

			//	double SM = 0.0;
			
			double Kap = MC.h1/MC.h2;

			double Teta = 0;
			double A =0.0;
			double B = 0.0;
			
			double Eps1out_tem = 0;
			double Eps1in_tem = 0;

			
			
			// coefficients for combined equations
			double A1 = 0.0;
			double B1 = 0.0;
			double V1 = 0.0;
			double V2 = 0.0;
			double W1 = 0.0;
			double W2 = 0.0;
			double K1 = 0.0;
			double K2 = 0.0; 
			double C1 = 0.0;
			


			// initial state:
			Tem = MC.Tmin;
			Fi = 1.0;
			M = 0.0;
			hM = 0.0;
			Sig1= 0.0;
			S1= 0.0;
			Sig2= 0.0;
			S2= 0.0;
			// ! hardening due to previous deformation must be taken into account!
			double SY1out_new=MC.SY1, SY1in_new =MC.SY1, SY2in_new =MC.SY2, SY2out_new=MC.SY2; //yield stress, 
			double SY1Ain_new= MC.SY1A, SY1Aout_new= MC.SY1A;
			// ! hardening due to previous deformation must be taken into account!

		
			Eps1out=0.0;
			Eps1in=0.0;
			Eps2in=0.0;
			Eps2out=0.0;
			Eps1outF=0.0;
			Eps1inF=0.0;
			Eps1outMP=0.0;
			Eps1inMP=0.0;
			EpsElastic=0.0;

						
			double Sig1out = 0.0;
			double Sig1in = 0.0;
			double Sig2in = 0.0;
			double Sig2out = 0.0;
			double CheckPlaneCond = 0.0;
			double Check_h_PlaneCond = 0.0;

			double DetGlav = 0.0;
			double DetX = 0.0;
			double DetY = 0.0;

			Hev1out = 0.0; //Hev(|Sig1+S1|)
			Hev1in = 0.0; //Hev(|Sig1-S1|)
			Hev2in = 0.0; //Hev(|Sig2+S2|)
			Hev2out = 0.0; //Hev(|Sig2-S2|

			Hev1out_old = 0.0; //Hev(|Sig1+S1|)
			Hev1in_old = 0.0; //Hev(|Sig1-S1|)
			Hev2in_old = 0.0; //Hev(|Sig2+S2|)
			Hev2out_old = 0.0; //Hev(|Sig2-S2|

			Hev1out_new = 0.0; //Hev(|Sig1+S1|)
			Hev1in_new = 0.0; //Hev(|Sig1-S1|)
			Hev2in_new = 0.0; //Hev(|Sig2+S2|)
			Hev2out_new = 0.0; //Hev(|Sig2-S2|)

			double Sig1_out_Yield;
			double Sig1_in_Yield;
			double Sig2_in_Yield;
			double Sig2_out_Yield;


			CheckPlaneCond = 0.0; // parameter for contol of the plane cross-section hypothesis

			res = fopen("res.dat","w");
			res_sm = fopen("res_sm.dat","w");
			fprintf(res,"# Tem\tM\tSig1\tS1\tSig2\tS2\tSig1out\tSig1in\tSig2in\tSig2out\tEps1out\tEps1in\tEps2in\tEps2out\tFi\tSigElastic\tEpsElastic\tCheckPlain");
			fprintf(res, "\n%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg", 
							Tem, M, Sig1*1e-6, S1*1e-6, Sig2*1e-6,S2*1e-6,Sig1out*1e-6,Sig1in*1e-6,Sig2in*1e-6,Sig2out*1e-6,
							Eps1out*1e2,Eps1in*1e2,Eps2in*1e2,Eps2out*1e2,Fi,6.0*M/MC.h/MC.h*1e-6,EpsElastic*1e2, MC.h2*(Eps1out-Eps1in)+MC.h1*(Eps2out-Eps2in), Eps1out_tem*1e2, Eps1in_tem*1e2);
			
										
//__________________________________loading___________________________________
				//hM = 1.0;


					Sig1_out_Yield = MC.SY1;
				Sig1_in_Yield = MC.SY1;
				Sig2_in_Yield = MC.SY2;
				Sig2_out_Yield = MC.SY2;

				// СПФ РАСТЯГИВАЕТСЯ + ЕСТЬ НАЧАЛЬНАЯ КРИВИЗНА
				//при растяжении наружного волокна кристалличсекого слоя возникают напряжения Sig1out=Sig1+S1.Ниже вычисляются слагаемые Sig1 и S1
				//Sig1 = Sig1_out_Yield + ((2*MC.r_l-(MC.l_n/MC.rad)*MC.h1)/((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l)-2*Sig1_out_Yield*D1M_L)*MC.H1/2 ;
				Sig1 = 0.5*(Sig1_out_Yield + (MC.r_l/((MC.l_n/MC.rad)*(MC.h1+MC.h2 + MC.rad)-MC.r_l)-Sig1_out_Yield*D1M_L)*MC.H1);
				//S1 = ((MC.l_n/MC.rad)*MC.h1/((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l))*MC.H1/2 ;
				S1 = 0.5*(Sig1_out_Yield + (MC.r_l/((MC.l_n/MC.rad)*(MC.h1+MC.h2 + MC.rad)-MC.r_l)-Sig1_out_Yield*D1M_L)*MC.H1);
				Sig2 = 0;
				S2 = 0; 	

				Sig1out = Sig1+S1;
				Sig1in = Sig1-S1;
				Sig2in = Sig2+S2;
				Sig2out = Sig2-S2;

				M = MC.h1*MC.h1*((3*Sig1+S1)/6); // изгибающий момент
				// деформация, соответствующая данному напряжению
				Eps1out = MC.r_l/((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l);
				//Eps1in = (MC.r_l-(MC.l_n/MC.rad)*MC.h1)/((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l);
				Eps1in = 0.0;
				Eps2in = 0;
				Eps2out =0;

				//*************************
				double Eps1outMax = Eps1out;
				double Eps1inMax = Eps1in;
				//*************************

				Hev1out = Hev(fabs(Sig1out)-MC.SY1); //Hev(|Sig1+S1|)			
				Hev1in = Hev(fabs(Sig1in)-MC.SY1); //Hev(|Sig1-S1|)				
				Hev2in = Hev(fabs(Sig2in)-MC.SY2); //Hev(|Sig2+S2|)				
				Hev2out = Hev(fabs(Sig2out)-MC.SY2); //Hev(|Sig2-S2|)

				if (Sig1out <= 0){ 
							Sig1_out_Yield = (-1.0)*MC.SY1;}
						if (Sig1in <= 0){ 
							Sig1_in_Yield = (-1.0)*MC.SY1;}
						if (Sig2in <= 0){ 
							Sig2_in_Yield = (-1.0)*MC.SY2;}
						if (Sig2out <= 0){ 
							Sig2_out_Yield = (-1.0)*MC.SY2;}

						//Sig1 = Sig1_out_Yield + ((2*MC.r_l-(MC.l_n/MC.rad)*MC.h1)/((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l)-2*Sig1_out_Yield*D1M_L)*MC.H1/2 ;
						Sig1 = 0.5*(Sig1_out_Yield + (MC.r_l/((MC.l_n/MC.rad)*(MC.h1+MC.h2 + MC.rad)-MC.r_l)-Sig1_out_Yield*D1M_L)*MC.H1);
						//S1 = ((MC.l_n/MC.rad)*MC.h1/((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l))*MC.H1/2 ;
						S1 = 0.5*(Sig1_out_Yield + (MC.r_l/((MC.l_n/MC.rad)*(MC.h1+MC.h2 + MC.rad)-MC.r_l)-Sig1_out_Yield*D1M_L)*MC.H1);
						Sig2 = 0;
						S2 = 0; 

						Sig1out = Sig1+S1;
						Sig1in = Sig1-S1;
						Sig2in = Sig2+S2;
						Sig2out = Sig2-S2;

						M = MC.h1*MC.h1*((3*Sig1+S1)/6);

						Eps1out = MC.r_l/((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l);
						//Eps1in = (MC.r_l-(MC.l_n/MC.rad)*MC.h1)/((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l);
						Eps1in = 0;
						Eps2in = 0;
						Eps2out =0;

						 

						EpsElastic = ((MC.l_n/MC.rad)*(MC.h1+MC.h2-MC.rad*Eps2out)-MC.r_l+Eps1out*((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l))/((MC.l_n/MC.rad)*(2*MC.rad+MC.h1+MC.h2+Eps2out*MC.rad)-MC.r_l+Eps1out*((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l)) ;

						fprintf(res, "\n%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg", 
							Tem, M, Sig1*1e-6, S1*1e-6, Sig2*1e-6,S2*1e-6,Sig1out*1e-6,Sig1in*1e-6,Sig2in*1e-6,Sig2out*1e-6,
							Eps1out*1e2,Eps1in*1e2,Eps2in*1e2,Eps2out*1e2,Fi,6.0*M/MC.h/MC.h*1e-6,EpsElastic*1e2, MC.h2*(Eps1out-Eps1in)+MC.h1*(Eps2out-Eps2in), Eps1out_tem*1e2,Eps1in_tem*1e2);
//___________________________unloading____________________________________________
		
				Hev1out_old = 0.0;
				Hev1in_old = 0.0;
				Hev2in_old = 0.0;
				Hev2out_old = 0.0;


				hM = -0.001;//-1e-3;//-0.01;

			while(M>0.0){
				double M_new=M+hM; if(M_new<1e-12)M_new=0.0;
				bool New_Old_Compar = false; // flag for cheking of conditions of plasticity: 
	
				A1 = ((Podat1U*(Hev1out_old-Hev1in_old) - Kap*Kap*Podat2*(Hev2out_old-Hev2in_old)))/
						Kap/(2.0*D2 + Podat2*(Hev2out_old+Hev2in_old));
				B1 = (2.0/MC.E1M_U + Podat1U*(Hev1out_old+Hev1in_old))/Kap/(2.0*D2 + Podat2*(Hev2out_old+Hev2in_old)); 
				V1 = 3.0*Kap*(Kap+1.0)+A1;
				W1 = Kap*Kap+B1;
				K1 = 6.0*hM/MC.h2/MC.h2;

				V2 = 1.0/MC.E1M_U + Podat1U*Hev1in_old + (D2+Podat2*Hev2in_old)*(Kap-A1);
				W2 = -B1*(D2+Podat2*Hev2in_old) - 1.0/MC.E1M_U - Podat1U*Hev1in_old; 
				K2 = 0.0;

		//		Kramer(V1,V2,W1,W2,K1,K2, hSig1,hS1);


				DetGlav = V1*W2 - W1*V2;
				DetX = K1*W2 - K2*W1;
				DetY = V1*K2 - V2*K1;


				hSig1 = DetX/DetGlav;
				hS1 = DetY/DetGlav;

				hSig2 = -Kap*hSig1;
				hS2 = A1*hSig1 + B1*hS1;

/*
				Sig1out_new = Sig1+hSig1+S1+hS1;
				Sig1in_new = Sig1+hSig1-S1-hS1;
				Sig2in_new = Sig2+hSig2+S2+hS2;
				Sig2out_new = Sig2+hSig2-S2-hS2;
				
				Hev1out_new = Hev(fabs(Sig1out_new)-SY1out_new)*Hev(fabs(Sig1out_new)-fabs(Sig1out)); //Hev(|Sig1+S1|)
				Hev1in_new = Hev(fabs(Sig1in_new)-SY1in_new)*Hev(fabs(Sig1in_new)-fabs(Sig1in)); //Hev(|Sig1-S1|)
				Hev2in_new = Hev(fabs(Sig2in_new)-SY2in_new)*Hev(fabs(Sig2in_new)-fabs(Sig2in)); //Hev(|Sig2+S2|)
				Hev2out_new = Hev(fabs(Sig2out_new)-SY2out_new)*Hev(fabs(Sig2out_new)-fabs(Sig2out)); //Hev(|Sig2-S2|)
			
	//				if (Hev1out_new != Hev1out_old) {Hev1out_old = Hev1out_new; New_Old_Compar = true;}
	//				if (Hev1in_new  != Hev1in_old)  {Hev1in_old  = Hev1in_new;  New_Old_Compar = true;}
	//				if (Hev2out_new != Hev2out_old) {Hev2out_old = Hev2out_new; New_Old_Compar = true;}
	//				if (Hev2in_new  != Hev2in_old)  {Hev2in_old  = Hev2in_new;  New_Old_Compar = true;}
//
					if (New_Old_Compar){ //if something changed, we calculate stress again
						A1 = Podat1U*(Hev1out_old-Hev1in_old) - Kap*Kap*Podat2*(Hev2out_old-Hev2in_old)/
								Kap/(2.0*D2 + Podat2*(Hev2out_old+Hev2in_old));
						B1 = (2.0/MC.E1M_U + Podat1U*(Hev1out_old+Hev1in_old))/Kap/(2.0*D2 + Podat2*(Hev2out_old+Hev2in_old)); 
						V1 = 3.0*Kap*(Kap+1.0)+A1;
						W1 = Kap*Kap+B1;
						K1 = 6.0*hM/MC.h2/MC.h2;

						V2 = 1.0/MC.E1M_U + Podat1U*Hev1in_old + (D2+Podat2*Hev2in_old)*(Kap-A1);
						W2 = -B1*(D2+Podat2*Hev2in_old) - 1.0/MC.E1M_U - Podat1U*Hev1in_old; 
						K2 = 0.0;

				//		Kramer(V1,V2,W1,W2,K1,K2, hSig1,hS1);

						DetGlav = V1*W1 - W1*V2;
						DetX = K1*W2 - K2*W1;
						DetY = V1*K2 - V2*K1;

						hSig1 = DetX/DetGlav;
						hS1 = DetY/DetGlav;

						hSig2 = -Kap*hSig1;
						hS2 = A1*hSig1 + B1*hS1;
					}
				

				//A1 = (D1M_U+Podat1U*Hev1in_old)/(D2+Podat2*Hev2in_old) + Kap;
				//B1 = (-D1M_U-Podat1U*Hev1in_old)/(D2+Podat2*Hev2in_old); 
				//V1 = 3.0*Kap*(Kap+1.0)+A1;
				//W1 = Kap*Kap+B1;
				//K1 = 6.0*hM/MC.h2/MC.h2;

				//V2 = Podat1U*(Hev1out_old - Hev1in_old) - Kap*Kap*Podat2*(Hev2out_old - Hev2in_old) 
				//	+ Kap*A1*(-2.0*D2+Podat2*(Hev2out_old + Hev2in_old));
				//W2 = 2.0*D1M_U + Podat1U*(Hev1out_old + Hev1in_old) 
				//	 +Kap*B1*(-2.0*D2+Podat2*(Hev2out_old - Hev2in_old)); 
				//K2 = 0.0;

				//Kramer(V1,V2,W1,W2,K1,K2, hSig1,hS1);
				//hSig2 = -Kap*hSig1;
				//hS2 = A1*hSig1 + B1*hS1;

				//Sig1out_new = Sig1+hSig1+S1+hS1;
				//Sig1in_new = Sig1+hSig1-S1-hS1;
				//Sig2in_new = Sig2+hSig2+S2+hS2;
				//Sig2out_new = Sig2+hSig2-S2-hS2;
				//
				//Hev1out_new = Hev(fabs(Sig1out_new)-MC.SY1)*Hev(fabs(Sig1out_new)-fabs(Sig1out)); //Hev(|Sig1+S1|)
				//Hev1in_new = Hev(fabs(Sig1in_new)-MC.SY1)*Hev(fabs(Sig1in_new)-fabs(Sig1in)); //Hev(|Sig1-S1|)
				//Hev2in_new = Hev(fabs(Sig2in_new)-MC.SY2)*Hev(fabs(Sig2in_new)-fabs(Sig2in)); //Hev(|Sig2+S2|)
				//Hev2out_new = Hev(fabs(Sig2out_new)-MC.SY2)*Hev(fabs(Sig2out_new)-fabs(Sig2out)); //Hev(|Sig2-S2|)
			
				//	if (Hev1out_new != Hev1out_old) {Hev1out_old = Hev1out_new; New_Old_Compar = true;}
				//	if (Hev1in_new  != Hev1in_old)  {Hev1in_old  = Hev1in_new;  New_Old_Compar = true;}
				//	if (Hev2out_new != Hev2out_old) {Hev2out_old = Hev2out_new; New_Old_Compar = true;}
				//	if (Hev2in_new  != Hev2in_old)  {Hev2in_old  = Hev2in_new;  New_Old_Compar = true;}

				//	if (New_Old_Compar){ //if something has changed, we calculate stress again
				//		A1 = (D1M_U+Podat1U*Hev1in_old)/(D2+Podat2*Hev2in_old) + Kap;
				//		B1 = (-D1M_U-Podat1U*Hev1in_old)/(D2+Podat2*Hev2in_old); 
				//		V1 = 3.0*Kap*(Kap+1.0)+A1;
				//		W1 = Kap*Kap+B1;
				//		K1 = 6.0*hM/MC.h2/MC.h2;

				//		V2 = Podat1U*(Hev1out_old - Hev1in_old) - Kap*Kap*Podat2*(Hev2out_old - Hev2in_old) 
				//			+ Kap*A1*(-2.0*D2+Podat2*(Hev2out_old + Hev2in_old));
				//		W2 = 2.0*D1M_U + Podat1U*(Hev1out_old + Hev1in_old) 
				//			 +Kap*B1*(-2.0*D2+Podat2*(Hev2out_old + Hev2in_old)); 
				//		K2 = 0.0;

				//		Kramer(V1,V2,W1,W2,K1,K2, hSig1,hS1);
				//		hSig2 = -Kap*hSig1;
				//		hS2 = A1*hSig1 + B1*hS1;
				//	}				
*/
				Eps1out += (hSig1+hS1)/MC.E1M_U + Podat1L*(hSig1+hS1)*Hev1out_old;
				Eps1in += (hSig1-hS1)/MC.E1M_U + Podat1L*(hSig1-hS1)*Hev1in_old;
				Eps2in += (hSig2+hS2)/MC.E2 + Podat2*(hSig2+hS2)*Hev2in_old;
				Eps2out += (hSig2-hS2)/MC.E2 + Podat2*(hSig2-hS2)*Hev2out_old;
				
				Sig1+=hSig1; // for output in the file
				S1+=hS1;
				Sig2+=hSig2;
				S2+=hS2;

				Sig1out = Sig1+S1; // stresses in layers
				Sig1in = Sig1-S1;
				Sig2in = Sig2+S2;
				Sig2out = Sig2-S2;

				M = M_new;

				

				EpsElastic = ((MC.l_n/MC.rad)*(MC.h1+MC.h2-MC.rad*Eps2out)-MC.r_l+Eps1out*((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l))/((MC.l_n/MC.rad)*(2*MC.rad+MC.h1+MC.h2+Eps2out*MC.rad)-MC.r_l+Eps1out*((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l)) ;
				//EpsElastic = (Eps1out - Eps2out)/2.0; ЗАМЕНИЛИ УПУРУГУЮ ДЕФОРМАЦИЮ 10 МАЯ 2017 СМ В ТЕТРАДИ

			fprintf(res, "\n%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg", 
							Tem, M, Sig1*1e-6, S1*1e-6, Sig2*1e-6,S2*1e-6,Sig1out*1e-6,Sig1in*1e-6,Sig2in*1e-6,Sig2out*1e-6,
							Eps1out*1e2,Eps1in*1e2,Eps2in*1e2,Eps2out*1e2,Fi,6.0*M/MC.h/MC.h*1e-6,EpsElastic*1e2, MC.h2*(Eps1out-Eps1in)+MC.h1*(Eps2out-Eps2in),Eps1out_tem*1e2 , Eps1in_tem*1e2);
			}//while M>0.0
			
			double EpsUnload = EpsElastic;
			fprintf(res_sm, "#EpsT\th1\th2\th\tk_h\tEpsUnload\tEps1out_heat\tTWSM\tKrBimet\n");
			fprintf(res_sm, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg", MC.EpsT*1e2, MC.h1, MC.h2, MC.h, MC.k_h, EpsUnload*1e2);


			//	if (fabs(Sig1out)>SY1out_new) SY1out_new = fabs(Sig1out); // new yield stress in out-layer of TiNi after unloading
			//	if (fabs(Sig1in)>SY1in_new) SY1in_new = fabs(Sig1in); // new yield stress in in-layer of TiNi after unloading
				SY1out_new = MC.SY1A; // new yield stress in out-layer of TiNi after unloading
				SY1in_new = MC.SY1A; // new yield stress in in-layer of TiNi after unloading

				if (fabs(Sig2in)>SY2in_new) SY2in_new = fabs(Sig2in); // new yield stress in in-layer of steel after unloading
				if (fabs(Sig2out)>SY2out_new) SY2out_new = fabs(Sig2out); // new yield stress in out-layer of steel after unloading

				//Eps1outF = Eps1out; // phase deformation in out-layer
				//Eps1inF = Eps1in; // phase deformation in in-layer
				Eps1outF = (Eps1outMax - Sig1_out_Yield/MC.E1M_U);
				Eps1inF = (Eps1inMax - Sig1_in_Yield/MC.E1M_U);
				//Eps1outF = Eps1out - MC.SY1/MC.E1M_U;
				//Eps1inF = Eps1in - MC.SY1/MC.E1M_U;
				Eps1outMP= (1-MC.Kr)*Eps1out; // irreversible deformation after heating  in out-layer
				Eps1inMP=(1-MC.Kr)*Eps1in; // irreversible deformation after heating  in in-layer

				



//_________________________heating-SM________________________________________________________________
		
		hTem = 0.0005;
		double Alf1_old = MC.Alf1A*(1.0-Fi)+Fi*MC.Alf1M;
		double D1_old = 1.0/(MC.E1A*(1.0-Fi)+Fi*MC.E1M_L);
		double Eps0out = MC.Kr*Eps1outF; // coefficient for calc. of phase deformation in out-layer
		double Eps0in = MC.Kr*Eps1inF; // coefficient for calc. of phase deformation in in-layer
		


	//Hevisid's functions for martensite and austenite	
		double Hev1outM_old = 0.0; //Hev(|Sig1+S1|-SY1)
		double Hev1outA_old = 0.0; //Hev(|Sig1+S1|)-SY1A)
		double Hev1inM_old = 0.0; //Hev(|Sig1-S1|-SY1)
		double Hev1inA_old = 0.0; //Hev(|Sig1-S1|)-SY1A)
		double Hev1outM_new = 0.0; //Hev(|Sig1+S1|-SY1)
		double Hev1outA_new = 0.0; //Hev(|Sig1+S1|)-SY1A)
		double Hev1inM_new = 0.0; //Hev(|Sig1-S1|-SY1)
		double Hev1inA_new = 0.0; //Hev(|Sig1-S1|)-SY1A)

		while (Tem < MC.Tmax){
			double Tem_new=Tem+hTem;
			bool New_Old_Compar = false; // flag for cheking of conditions of plasticity:

			Teta = Tem_new - MC.Tmin;
			  //Teta = Tem - MC.Ms;
			Pt = Tem + Fi*(MC.Af-MC.As) - MC.Af;
			hFi = -hTem/(MC.Af - MC.As)*Hev(Fi,0)*Hev(Pt,1);
			//double Fi_new = Fi+hFi; // переменная для сравнения фазы 
			if(Fi+hFi<0.0) hFi=-Fi;
			
		//		D1_new = 1.0/(MC.E1A*(1.0-Fi)+Fi*MC.E1M_L);
	//		hD1 = D1_new-D1_old;
			
			double Alf1_new = MC.Alf1A*(1.0-Fi)+Fi*MC.Alf1M;
			double hAlf1 = Alf1_new - Alf1_old;

			V1 = Fi*(D1M_L+Podat1L*Hev1inM_old)+(1.0-Fi)*(D1A+Podat1A*Hev1inA_old)				
				+ Kap*(3.0*Kap+4.0)*(D2 + Podat2*Hev2in_old);
			
			W1 = -Fi*(D1M_L+Podat1L*Hev1inM_old)-(1.0-Fi)*(D1A+Podat1A*Hev1inA_old)
				+ Kap*Kap*(D2 + Podat2*Hev2in_old);
			
			V2 = Fi*Podat1L*(Hev1outM_old-Hev1inM_old) + (1.0-Fi)*Podat1A*(Hev1outA_old-Hev1inA_old)
				- Kap*Kap*Podat2*(Hev2out_old-Hev2in_old)
				 +3.0*Kap*Kap*(Kap+1)*(2.0*D2+Podat2*(Hev2out_old+Hev2in_old));				   
			
			W2 = Fi*(2.0*D1M_L+Podat1L*(Hev1outM_old+Hev1inM_old)) + (1.0-Fi)*(2.0*D1A+Podat1A*(Hev1outA_old+Hev1inA_old)) 
				 + Kap*Kap*Kap*(2.0*D2 + Podat2*(Hev2out_old+Hev2in_old));
			
			K1 = hFi*((-D1M_L-Podat1L*Hev1inM_old+D1A+Podat1A*Hev1inA_old)*Sig1in - Eps0in)
				-hTem*(MC.Alf1M-MC.Alf1A)*hFi
				/*- (Alf1_old - MC.Alf2)*hTem - Tem*hAlf1*/; // Sig1in = Sig1-S1 //+MC.Alf2*hTem-Alf1_new*hTem
			
			K2 = -hFi*(2.0*S1*(D1M_L-D1A) + Sig1out*(Podat1L*Hev1outM_old-Podat1A*Hev1outA_old) 
				 - Sig1in*(Podat1L*Hev1inM_old-Podat1A*Hev1inA_old) + Eps0out - Eps0in);

			Kramer(V1,V2,W1,W2,K1,K2, hSig1,hS1);
			hSig2 = -Kap*hSig1;
			hS2 = -3.0*Kap*(Kap+1.0)*hSig1 - Kap*Kap*hS1;
			
			Sig1out_new = Sig1+hSig1+S1+hS1;
			Sig1in_new = Sig1+hSig1-S1-hS1;
			Sig2in_new = Sig2+hSig2+S2+hS2;
			Sig2out_new = Sig2+hSig2-S2-hS2;
			
/*
	// now we try to account the plastic deformation of TiNi at heating (only for austenite)
			Hev1outM_new = 0.0;//Hev(fabs(Sig1out_new)-SY1out_new)*Hev(fabs(Sig1out_new)-fabs(Sig1out)); //Hev(|Sig1+S1|)
			Hev1inM_new = 0.0;//Hev(fabs(Sig1in_new)-SY1in_new)*Hev(fabs(Sig1in_new)-fabs(Sig1in)); //Hev(|Sig1-S1|)
//			Hev1outA_new = Hev(fabs(Sig1out_new)-MC.SY1A)*Hev(fabs(Sig1out_new)-fabs(Sig1out)); //Hev(|Sig1+S1|)
//			Hev1inA_new = Hev(fabs(Sig1in_new)-MC.SY1A)*Hev(fabs(Sig1in_new)-fabs(Sig1in)); //Hev(|Sig1-S1|)
	
			Hev2in_new = Hev(fabs(Sig2in_new)-SY2in_new)*Hev(fabs(Sig2in_new)-fabs(Sig2in)); //Hev(|Sig2+S2|)
			Hev2out_new = Hev(fabs(Sig2out_new)-SY2out_new)*Hev(fabs(Sig2out_new)-fabs(Sig2out)); //Hev(|Sig2-S2|)

			if (Hev1outM_new != Hev1outM_old) {Hev1outM_old = Hev1outM_new; New_Old_Compar = true;}
			if (Hev1inM_new  != Hev1inM_old)  {Hev1inM_old  = Hev1inM_new;  New_Old_Compar = true;}
//			if (Hev1outA_new != Hev1outA_old) {Hev1outA_old = Hev1outA_new; New_Old_Compar = true;}
//			if (Hev1inA_new  != Hev1inA_old)  {Hev1inA_old  = Hev1inA_new;  New_Old_Compar = true;}			
			if (Hev2out_new != Hev2out_old) {Hev2out_old = Hev2out_new; New_Old_Compar = true;}
			if (Hev2in_new  != Hev2in_old)  {Hev2in_old  = Hev2in_new;  New_Old_Compar = true;}

			if (New_Old_Compar){ //if something changed, we calculate stress again

			V1 = Fi*(D1M_L+Podat1L*Hev1inM_old)+(1.0-Fi)*(D1A+Podat1A*Hev1inA_old)				
				+ Kap*(3.0*Kap+4.0)*(D2 + Podat2*Hev2in_old);
			
			W1 = -Fi*(D1M_L+Podat1L*Hev1inM_old)-(1.0-Fi)*(D1A+Podat1A*Hev1inA_old)
				+ Kap*Kap*(D2 + Podat2*Hev2in_old);
			
			V2 = Fi*Podat1L*(Hev1outM_old-Hev1inM_old) + (1.0-Fi)*Podat1A*(Hev1outA_old-Hev1inA_old)
				- Kap*Kap*Podat2*(Hev2out_old-Hev2in_old)
				 +3.0*Kap*Kap*(Kap+1)*(2.0*D2+Podat2*(Hev2out_old+Hev2in_old));				   
			
			W2 = Fi*(2.0*D1M_L+Podat1L*(Hev1outM_old+Hev1inM_old)) + (1.0-Fi)*(2.0*D1A+Podat1A*(Hev1outA_old+Hev1inA_old)) 
				 + Kap*Kap*Kap*(2.0*D2 + Podat2*(Hev2out_old+Hev2in_old));
			
			K1 = hFi*((-D1M_L-Podat1L*Hev1inM_old+D1A+Podat1A*Hev1inA_old)*Sig1in - Eps0in)
				- (Alf1_old - MC.Alf2)*hTem - Tem*hAlf1; // Sig1in = Sig1-S1
			
			K2 = -hFi*(2.0*S1*(D1M_L-D1A) + Sig1out*(Podat1L*Hev1outM_old-Podat1A*Hev1outA_old) 
				 - Sig1in*(Podat1L*Hev1inM_old-Podat1A*Hev1inA_old) - Eps0out + Eps0in);

			Kramer(V1,V2,W1,W2,K1,K2, hSig1,hS1);
			hSig2 = -Kap*hSig1;
			hS2 = -3.0*Kap*(Kap+1.0)*hSig1 - Kap*Kap*hS1;
			}
*/
		//Update_step_heating1:
			Eps1out_tem =  hTem*(MC.Alf1M-MC.Alf1A)*hFi; //Alf1_new*hTem +
			Eps1in_tem = hTem*(MC.Alf1M-MC.Alf1A)*hFi; //Alf1_new*hTem + 
			A = hFi*Sig1out*(D1M_L + Podat1L*Hev1outM_new - D1A - Podat1A*Hev1outA_new) + hFi*Eps0out
						+ (hSig1+hS1)*(Fi*(D1M_L + Podat1L*Hev1outM_new) + (1.0-Fi)*(D1A + Podat1A*Hev1outA_new));
			B = hFi*Sig1in*(D1M_L + Podat1L*Hev1inM_new - D1A - Podat1A*Hev1inA_new) + hFi*Eps0in
						+ (hSig1-hS1)*(Fi*(D1M_L + Podat1L*Hev1inM_new) + (1.0-Fi)*(D1A + Podat1A*Hev1inA_new));
					
			Eps1out += /*hFi*Sig1out*(D1M_L + Podat1L*Hev1outM_new - D1A - Podat1A*Hev1outA_new) + hFi*Eps0out
						+ (hSig1+hS1)*(Fi*(D1M_L + Podat1L*Hev1outM_new) + (1.0-Fi)*(D1A + Podat1A*Hev1outA_new))*/ A + Eps1out_tem;			
						//+Alf1_new*hTem + Teta*(MC.Alf1M-MC.Alf1A)*hFi/*+ Alf1_old*hTem + Tem*hAlf1*/;
			Eps1in += /*hFi*Sig1in*(D1M_L + Podat1L*Hev1inM_new - D1A - Podat1A*Hev1inA_new) + hFi*Eps0in
						+ (hSig1-hS1)*(Fi*(D1M_L + Podat1L*Hev1inM_new) + (1.0-Fi)*(D1A + Podat1A*Hev1inA_new))*/ B	+Eps1in_tem;
						//+Alf1_new*hTem + Teta*(MC.Alf1M-MC.Alf1A)*hFi
						/*+ Alf1_old*hTem + Tem*hAlf1*/;
			Eps2in += (hSig2+hS2)/MC.E2 + Podat2*(hSig2+hS2)*Hev2in_new + MC.Alf2*hTem;
			Eps2out += (hSig2-hS2)/MC.E2 + Podat2*(hSig2-hS2)*Hev2out_new + MC.Alf2*hTem;


			Eps1outF += Eps0out*hFi;
			Eps1inF += Eps0in*hFi;

			//Eps1out_tem = Alf1_new*hTem + Teta*(MC.Alf1M-MC.Alf1A)*hFi;
			//Eps1in_tem = Alf1_new*hTem + Teta*(MC.Alf1M-MC.Alf1A)*hFi;

			Sig1+=hSig1;
			S1+=hS1;
			Sig2+=hSig2;
			S2+=hS2;
			Sig1out = Sig1+S1;
			Sig1in = Sig1-S1;
			Sig2in = Sig2+S2;
			Sig2out = Sig2-S2;

			if((Fi>0)&&((Fi+hFi)==0)){
				Eps0outAf = Eps1out;
				Eps0inAf = Eps1in;
				//std::cout<<"Eps0outAf=" <<Eps0outAf<<std::endl;
			}

			Fi += hFi;
			Tem = Tem_new;
	//		D1_old = D1_new;
			Alf1_old = Alf1_new;
			//EpsElastic = (Eps1out - Eps2out)/2.0;
			EpsElastic = ((MC.l_n/MC.rad)*(MC.h1+MC.h2-MC.rad*Eps2out)-MC.r_l+Eps1out*((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l))/((MC.l_n/MC.rad)*(2*MC.rad+MC.h1+MC.h2+Eps2out*MC.rad)-MC.r_l+Eps1out*((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l)) ;
	//		if (Tem == MC.Af) SM = Eps1out;

			fprintf(res, "\n%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg", 
							Tem, M, Sig1*1e-6, S1*1e-6, Sig2*1e-6,S2*1e-6,Sig1out*1e-6,Sig1in*1e-6,Sig2in*1e-6,Sig2out*1e-6,
							Eps1out*1e2,Eps1in*1e2,Eps2in*1e2,Eps2out*1e2,Fi,6.0*M/MC.h/MC.h*1e-6,EpsElastic*1e2, MC.h2*(Eps1out-Eps1in)+MC.h1*(Eps2out-Eps2in),Eps1out_tem*1e2, Eps1in_tem*1e2 );
		}//while (Tem < MC.Tmax)

			double EpsHeating = EpsElastic;
			double KrBimet = (EpsUnload - EpsHeating)/EpsUnload;

			fprintf(res_sm, "\t%lg", Eps1out*1e2);
			
			if (fabs(Sig1out)>MC.SY1A) SY1Aout_new = fabs(Sig1out); // new yield stress in out-layer of TiNi after unloading
			if (fabs(Sig1in)>MC.SY1A) SY1Ain_new = fabs(Sig1in); // new yield stress in in-layer of TiNi after unloading
			if (fabs(Sig2in)>SY2in_new) SY2in_new = fabs(Sig2in); // new yield stress in steel in-layer after 1-st heating
			if (fabs(Sig2out)>SY2out_new) SY2out_new = fabs(Sig2out); // new yield stress in steel out-layer after 1-st heating

			Eps1outF = MC.Kr*Eps1out;
			Eps1inF = MC.Kr*Eps1in;

//__________________________________TWSM-cooling______________________________


		//hTem = (MC.Mf-MC.Ms)*MC.ETP/MC.E1A/5.0; //the step on temperarute must be less then (MC.Mf-MC.Ms)*MC.ETP/MC.E1A
		hTem = -0.0005;									   //to avoid instability;
		double Eps0out0 = MC.Lam1*Eps1out; // coefficient for calc. of phase deformation in out-layer
		double Eps0in0 = MC.Lam1*Eps1in; // coefficient for calc. of phase deformation in in-layer
	
		while (Tem > MC.Tmin){
			double Tem_new=Tem+hTem;
			bool New_Old_Compar = false; // flag for cheking of conditions of plasticity:

			Teta = Tem_new - MC.Tmin;

			//Teta = Tem - MC.Ms;

			Pt = MC.Ms-Fi*(MC.Ms-MC.Mf)-Tem;
			hFi = -hTem/(MC.Ms - MC.Mf)*Hev(1-Fi,0)*Hev(Pt,1);
			if(Fi+hFi>1.0) hFi=1.0-Fi;
			//Fi += hFi;
			double Fi_new = Fi+hFi;
			double D1_new = 1.0/(MC.E1A*(1.0-Fi_new)+Fi_new*MC.E1M_L);
            double hD1 = D1_new-D1_old;
			Eps0out = Eps0out0 + Sig1out/MC.ETP; // coefficient for calc. of phase deformation in out-layer
			Eps0in = Eps0in0 + Sig1in/MC.ETP; // coefficient for calc. of phase deformation in in-layer

			double Alf1_new = MC.Alf1A*(1.0-Fi_new)+Fi_new*MC.Alf1M;
			double hAlf1 = Alf1_new - Alf1_old;


			V1 = Fi_new*(D1M_L+Podat1L*Hev1inM_old)+(1.0-Fi_new)*(D1A+Podat1A*Hev1inA_old)				
				+ Kap*(3.0*Kap+4.0)*(D2 + Podat2*Hev2in_old);
			
			W1 = -Fi_new*(D1M_L+Podat1L*Hev1inM_old)-(1.0-Fi_new)*(D1A+Podat1A*Hev1inA_old)
				+ Kap*Kap*(D2 + Podat2*Hev2in_old);
			
			V2 = Fi_new*Podat1L*(Hev1outM_old-Hev1inM_old) + (1.0-Fi_new)*Podat1A*(Hev1outA_old-Hev1inA_old)
				- Kap*Kap*Podat2*(Hev2out_old-Hev2in_old)
				 +3.0*Kap*Kap*(Kap+1)*(2.0*D2+Podat2*(Hev2out_old+Hev2in_old));				   
			
			W2 = Fi_new*(2.0*D1M_L+Podat1L*(Hev1outM_old+Hev1inM_old)) + (1.0-Fi_new)*(2.0*D1A+Podat1A*(Hev1outA_old+Hev1inA_old)) 
				 + Kap*Kap*Kap*(2.0*D2 + Podat2*(Hev2out_old+Hev2in_old));
			
			K1 = hFi*((-D1M_L-Podat1L*Hev1inM_old+D1A+Podat1A*Hev1inA_old)*Sig1in - Eps0in)
				 - hTem*(MC.Alf1M-MC.Alf1A)*hFi   /*- (Alf1_old - MC.Alf2)*hTem - Tem*hAlf1*/; // Sig1in = Sig1-S1 //убрали + MC.Alf2 - Alf1_new*hTem //teta->hTem
			
			K2 = -hFi*(2.0*S1*(D1M_L-D1A) + Sig1out*(Podat1L*Hev1outM_old-Podat1A*Hev1outA_old) 
				 - Sig1in*(Podat1L*Hev1inM_old-Podat1A*Hev1inA_old) + Eps0out - Eps0in);

			Kramer(V1,V2,W1,W2,K1,K2, hSig1,hS1);
			hSig2 = -Kap*hSig1;
			hS2 = -3.0*Kap*(Kap+1.0)*hSig1 - Kap*Kap*hS1;


			Sig1out_new = Sig1+hSig1+S1+hS1;
			Sig1in_new = Sig1+hSig1-S1-hS1;
			Sig2in_new = Sig2+hSig2+S2+hS2;
			Sig2out_new = Sig2+hSig2-S2-hS2;
	/*		
		//	Hev1out_new = Hev(fabs(Sig1out_new)-SY1out_new)*Hev(fabs(Sig1out_new)-fabs(Sig1out)); //Hev(|Sig1+S1|)
		//	Hev1in_new = Hev(fabs(Sig1in_new)-SY2out_new)*Hev(fabs(Sig1in_new)-fabs(Sig1in)); //Hev(|Sig1-S1|)
		// there is no plastic deformation of TiNi during thermocycles
			Hev2in_new = Hev(fabs(Sig2in_new)-SY2in_new)*Hev(fabs(Sig2in_new)-fabs(Sig2in)); //Hev(|Sig2+S2|)
			Hev2out_new = Hev(fabs(Sig2out_new)-SY2out_new)*Hev(fabs(Sig2out_new)-fabs(Sig2out)); //Hev(|Sig2-S2|)

	//		if (Hev1out_new != Hev1out_old) {Hev1out_old = Hev1out_new; New_Old_Compar = true;}
	//		if (Hev1in_new  != Hev1in_old)  {Hev1in_old  = Hev1in_new;  New_Old_Compar = true;}
			if (Hev2out_new != Hev2out_old) {Hev2out_old = Hev2out_new; New_Old_Compar = true;}
			if (Hev2in_new  != Hev2in_old)  {Hev2in_old  = Hev2in_new;  New_Old_Compar = true;}

			if (New_Old_Compar){ //if something changed, we calculate stress again

				V1 = D1_old + Kap*(3.0*Kap+4.0)*(D2 + Podat2*Hev2in_old);
				W1 = -D1_old + Kap*Kap*(D2 + Podat2*Hev2in_old);
				V2 = - Kap*Kap*Podat2*(Hev2out_old-Hev2in_old)
					 +3.0*Kap*Kap*(Kap+1)*(2.0*D2+Podat2*(Hev2out_old+Hev2in_old));				  
				W2 = 2.0*D1_old + Kap*Kap*Kap*(2.0*D2 + Podat2*(Hev2out_old+Hev2in_old));
				K1 = -Sig1in*hD1 - Eps0in*hFi - (Alf1_old - MC.Alf2)*hTem - Tem*hAlf1;
				K2 = -2.0*S1*hD1 - (Eps0out-Eps0in)*hFi;

				Kramer(V1,V2,W1,W2,K1,K2, hSig1,hS1);
				hSig2 = -Kap*hSig1;
				hS2 = -3.0*Kap*(Kap+1.0)*hSig1 - Kap*Kap*hS1;
			}

		//Update_step_heating1:
		
*/          Eps1out_tem =  hTem*(MC.Alf1M-MC.Alf1A)*hFi; //убрали Alf1_new*hTem +
			Eps1in_tem = hTem*(MC.Alf1M-MC.Alf1A)*hFi; //убрали Alf1_new*hTem + 
			A = hFi*Sig1out*(D1M_L + Podat1L*Hev1outM_new - D1A - Podat1A*Hev1outA_new) + hFi*Eps0out
						+ (hSig1+hS1)*(Fi_new*(D1M_L + Podat1L*Hev1outM_new) + (1.0-Fi_new)*(D1A + Podat1A*Hev1outA_new));
			B = hFi*Sig1in*(D1M_L + Podat1L*Hev1inM_new - D1A - Podat1A*Hev1inA_new) + hFi*Eps0in
						+ (hSig1-hS1)*(Fi_new*(D1M_L + Podat1L*Hev1inM_new) + (1.0-Fi_new)*(D1A + Podat1A*Hev1inA_new));
			double C=A+Eps1out_tem;
			double Z = (hSig1+hS1)*Fi_new*D1M_L;
			Eps1out += /*hFi*Sig1out*(D1M_L + Podat1L*Hev1outM_new - D1A - Podat1A*Hev1outA_new) + hFi*Eps0out
						+ (hSig1+hS1)*(Fi*(D1M_L + Podat1L*Hev1outM_new) + (1.0-Fi)*(D1A + Podat1A*Hev1outA_new))*/ A +	Eps1out_tem;	
						//+ Alf1_new*hTem + Teta*(MC.Alf1M - MC.Alf1A)*hFi/*+ Alf1_old*hTem + Tem*hAlf1*/;
			Eps1in += /* hFi*Sig1in*(D1M_L + Podat1L*Hev1inM_new - D1A - Podat1A*Hev1inA_new) + hFi*Eps0in
						+ (hSig1-hS1)*(Fi*(D1M_L + Podat1L*Hev1inM_new) + (1.0-Fi)*(D1A + Podat1A*Hev1inA_new))*/ B	+ Eps1in_tem;	
						//+ Alf1_new*hTem + Teta*(MC.Alf1M - MC.Alf1A)*hFi/*+ Alf1_old*hTem + Tem*hAlf1*/;
			Eps2in += (hSig2+hS2)/MC.E2 + Podat2*(hSig2+hS2)*Hev2in_new + MC.Alf2*hTem;
			Eps2out += (hSig2-hS2)/MC.E2 + Podat2*(hSig2-hS2)*Hev2out_new + MC.Alf2*hTem;

				
			Eps1outF += Eps0out*hFi;
			Eps1inF += Eps0in*hFi;

			//Eps1out_tem = Alf1_new*hTem + Teta*(MC.Alf1M-MC.Alf1A)*hFi;
			//Eps1in_tem = Alf1_new*hTem + Teta*(MC.Alf1M-MC.Alf1A)*hFi;


			Sig1+=hSig1;
			S1+=hS1;
			Sig2+=hSig2;
			S2+=hS2;
			Sig1out = Sig1+S1;
			Sig1in = Sig1-S1;
			Sig2in = Sig2+S2;
			Sig2out = Sig2-S2;

			Tem = Tem_new;
			D1_old = D1_new;
			Alf1_old = Alf1_new;
			//EpsElastic = (Eps1out - Eps2out)/2.0;
			EpsElastic = ((MC.l_n/MC.rad)*(MC.h1+MC.h2-MC.rad*Eps2out)-MC.r_l+Eps1out*((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l))/((MC.l_n/MC.rad)*(2*MC.rad+MC.h1+MC.h2+Eps2out*MC.rad)-MC.r_l+Eps1out*((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l)) ;
			if ((Fi==0)&&(Fi_new>0)/*((Fi>0)&&(Fi_new-0)<=0.2)*//*Tem>=MC.Ms*/) {
		//		TWSM = EpsElastic;//Eps1out;//
				Eps1Fout_MS = Eps0outAf;//Eps1out;
				Eps1Fin_MS = Eps0inAf;//Eps1in;
			}	
			if ((Fi<1)&&(Fi_new>=1)/*(((Fi_new-1.0)<=0.2)&&((Fi+1.0)>0))*//*Tem <= MC.Mf*/) {
				Eps1Fout_MF = Eps1out;
				Eps1Fin_MF = Eps1in;
				TWSM = EpsElastic;
			}
			
			Fi = Fi_new;

			fprintf(res, "\n%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg", 
							Tem, M, Sig1*1e-6, S1*1e-6, Sig2*1e-6,S2*1e-6,Sig1out*1e-6,Sig1in*1e-6,Sig2in*1e-6,Sig2out*1e-6,
							Eps1out*1e2,Eps1in*1e2,Eps2in*1e2,Eps2out*1e2,Fi,6.0*M/MC.h/MC.h*1e-6,EpsElastic*1e2, MC.h2*(Eps1out-Eps1in)+MC.h1*(Eps2out-Eps2in),Eps1out_tem*1e2, Eps1in_tem*1e2);
//			TWSM = EpsElastic-TWSM;//Eps1out-TWSM;//only for Alf1=0 and Alf2=0!
//			fprintf(res_sm, "\t%lg\t%lg\n", TWSM*1e2, KrBimet);

}			if (fabs(Sig2in)>SY2in_new) SY2in_new = fabs(Sig2in); // new yield stress in steel in-layer after 1-st cooling
			if (fabs(Sig2out)>SY2out_new) SY2out_new = fabs(Sig2out); // new yield stress in steel out-layer after 1-st cooling

//__________________________________TWSM-heating______________________________



		//hTem = -hTem;
	      hTem = 0.0005;
		Eps0out0 = Eps1Fout_MF-Eps1Fout_MS; // coefficient for calc. of phase deformation in out-layer
		Eps0in0 = Eps1Fin_MF-Eps1Fin_MS; // coefficient for calc. of phase deformation in in-layer

			
		//	printf ("\nEps0out0=%lg, Eps0in0=%lg" 
		//		      ,Eps0out0,     Eps0in0);
		//	pause();
			


	//		Eps0out = MC.Lam1*Eps1out+MC.Lam2*Sig1out/MC.E1A; // coefficient for calc. of phase deformation in out-layer
	//		Eps0in = MC.Lam1*Eps1in+MC.Lam2*Sig1in/MC.E1A; // coefficient for calc. of phase deformation in in-layer
		while (Tem < MC.Tmax){
			double Tem_new=Tem+hTem;
			bool New_Old_Compar = false; // flag for cheking of conditions of plasticity:

			Teta = Tem_new - MC.Tmin; // не было ранее

			Pt = Tem + Fi*(MC.Af-MC.As) - MC.Af;
			hFi = -hTem/(MC.Af - MC.As)*Hev(Fi,0)*Hev(Pt,1);
			if(Fi+hFi<0.0) hFi=-Fi;
			Fi += hFi;
			double D1_new = 1.0/(MC.E1A*(1.0-Fi)+Fi*MC.E1M_L);
			double hD1 = D1_new-D1_old;
			double Alf1_new = MC.Alf1A*(1.0-Fi)+Fi*MC.Alf1M;
			double hAlf1 = Alf1_new - Alf1_old;
			Eps0out = Eps0out0;// coefficient for calc. of phase deformation in out-layer
			Eps0in = Eps0in0; // coefficient for calc. of phase deformation in in-layer

			V1 = Fi*(D1M_L+Podat1L*Hev1inM_old)+(1.0-Fi)*(D1A+Podat1A*Hev1inA_old)				
				+ Kap*(3.0*Kap+4.0)*(D2 + Podat2*Hev2in_old);
			
			W1 = -Fi*(D1M_L+Podat1L*Hev1inM_old)-(1.0-Fi)*(D1A+Podat1A*Hev1inA_old)
				+ Kap*Kap*(D2 + Podat2*Hev2in_old);
			
			V2 = Fi*Podat1L*(Hev1outM_old-Hev1inM_old) + (1.0-Fi)*Podat1A*(Hev1outA_old-Hev1inA_old)
				- Kap*Kap*Podat2*(Hev2out_old-Hev2in_old)
				 +3.0*Kap*Kap*(Kap+1)*(2.0*D2+Podat2*(Hev2out_old+Hev2in_old));				   
			
			W2 = Fi*(2.0*D1M_L+Podat1L*(Hev1outM_old+Hev1inM_old)) + (1.0-Fi)*(2.0*D1A+Podat1A*(Hev1outA_old+Hev1inA_old)) 
				 + Kap*Kap*Kap*(2.0*D2 + Podat2*(Hev2out_old+Hev2in_old));
			
			K1 = hFi*((-D1M_L-Podat1L*Hev1inM_old+D1A+Podat1A*Hev1inA_old)*Sig1in - Eps0in)
				 - hTem*(MC.Alf1M-MC.Alf1A)*hFi/*- (Alf1_old - MC.Alf2)*hTem - Tem*hAlf1*/; // Sig1in = Sig1-S1 //убрали + MC.Alf2 - Alf1_new*hTem
			
			K2 = -hFi*(2.0*S1*(D1M_L-D1A) + Sig1out*(Podat1L*Hev1outM_old-Podat1A*Hev1outA_old) 
				 - Sig1in*(Podat1L*Hev1inM_old-Podat1A*Hev1inA_old) + Eps0out - Eps0in);

			Kramer(V1,V2,W1,W2,K1,K2, hSig1,hS1);
			hSig2 = -Kap*hSig1;
			hS2 = -3.0*Kap*(Kap+1.0)*hSig1 - Kap*Kap*hS1;

			Sig1out_new = Sig1+hSig1+S1+hS1;
			Sig1in_new = Sig1+hSig1-S1-hS1;
			Sig2in_new = Sig2+hSig2+S2+hS2;
			Sig2out_new = Sig2+hSig2-S2-hS2;
	/*		
		//	Hev1out_new = Hev(fabs(Sig1out_new)-MC.SY1)*Hev(fabs(Sig1out_new)-fabs(Sig1out)); //Hev(|Sig1+S1|)
		//	Hev1in_new = Hev(fabs(Sig1in_new)-MC.SY1)*Hev(fabs(Sig1in_new)-fabs(Sig1in)); //Hev(|Sig1-S1|)
		// there is no plastic deformation of TiNi during thermocycles
			Hev2in_new = Hev(fabs(Sig2in_new)-SY2in_new)*Hev(fabs(Sig2in_new)-fabs(Sig2in)); //Hev(|Sig2+S2|)
			Hev2out_new = Hev(fabs(Sig2out_new)-SY2out_new)*Hev(fabs(Sig2out_new)-fabs(Sig2out)); //Hev(|Sig2-S2|)

	//		if (Hev1out_new != Hev1out_old) {Hev1out_old = Hev1out_new; New_Old_Compar = true;}
	//		if (Hev1in_new  != Hev1in_old)  {Hev1in_old  = Hev1in_new;  New_Old_Compar = true;}
			if (Hev2out_new != Hev2out_old) {Hev2out_old = Hev2out_new; New_Old_Compar = true;}
			if (Hev2in_new  != Hev2in_old)  {Hev2in_old  = Hev2in_new;  New_Old_Compar = true;}

			if (New_Old_Compar){ //if something changed, we calculate stress again

				V1 = D1_old + Kap*(3.0*Kap+4.0)*(D2 + Podat2*Hev2in_old);
				W1 = -D1_old + Kap*Kap*(D2 + Podat2*Hev2in_old);
				V2 = - Kap*Kap*Podat2*(Hev2out_old-Hev2in_old)
					 +3.0*Kap*Kap*(Kap+1)*(2.0*D2+Podat2*(Hev2out_old+Hev2in_old));				   
				W2 = 2.0*D1_old + Kap*Kap*Kap*(2.0*D2 + Podat2*(Hev2out_old+Hev2in_old));
				K1 = -Sig1in*hD1 - Eps0in*hFi - (Alf1_old - MC.Alf2)*hTem - Tem*hAlf1;
				K2 = -2.0*S1*hD1 - (Eps0out-Eps0in)*hFi;

				Kramer(V1,V2,W1,W2,K1,K2, hSig1,hS1);
				hSig2 = -Kap*hSig1;
				hS2 = -3.0*Kap*(Kap+1.0)*hSig1 - Kap*Kap*hS1;
			}

		//Update_step_heating1:
*/
			Eps1out_tem = hTem*(MC.Alf1M-MC.Alf1A)*hFi; //Alf1_new*hTem + 
			Eps1in_tem = hTem*(MC.Alf1M-MC.Alf1A)*hFi; // Alf1_new*hTem + 
			A = hFi*Sig1out*(D1M_L + Podat1L*Hev1outM_new - D1A - Podat1A*Hev1outA_new) + hFi*Eps0out
						+ (hSig1+hS1)*(Fi*(D1M_L + Podat1L*Hev1outM_new) + (1.0-Fi)*(D1A + Podat1A*Hev1outA_new));
			B = hFi*Sig1in*(D1M_L + Podat1L*Hev1inM_new - D1A - Podat1A*Hev1inA_new) + hFi*Eps0in
						+ (hSig1-hS1)*(Fi*(D1M_L + Podat1L*Hev1inM_new) + (1.0-Fi)*(D1A + Podat1A*Hev1inA_new));
			double C = A + Eps1out_tem;
			double Z = (hSig1-hS1)*Fi*D1M_L;

			Eps1out += /*hFi*Sig1out*(D1M_L + Podat1L*Hev1outM_new - D1A - Podat1A*Hev1outA_new) + hFi*Eps0out
						+ (hSig1+hS1)*(Fi*(D1M_L + Podat1L*Hev1outM_new) + (1.0-Fi)*(D1A + Podat1A*Hev1outA_new))*/ A +	Eps1out_tem;	
						//+ Alf1_new*hTem + Teta*(MC.Alf1M - MC.Alf1A)*hFi/*+ Alf1_old*hTem + Tem*hAlf1*/;
			Eps1in += /*hFi*Sig1in*(D1M_L + Podat1L*Hev1inM_new - D1A - Podat1A*Hev1inA_new) + hFi*Eps0in
						+ (hSig1-hS1)*(Fi*(D1M_L + Podat1L*Hev1inM_new) + (1.0-Fi)*(D1A + Podat1A*Hev1inA_new))*/ B	+Eps1in_tem;
						//+ Alf1_new*hTem + Teta*(MC.Alf1M - MC.Alf1A)*hFi/*+ Alf1_old*hTem + Tem*hAlf1*/;
			Eps2in += (hSig2+hS2)/MC.E2 + Podat2*(hSig2+hS2)*Hev2in_new + MC.Alf2*hTem;
			Eps2out += (hSig2-hS2)/MC.E2 + Podat2*(hSig2-hS2)*Hev2out_new + MC.Alf2*hTem;


			Eps1outF += Eps0out*hFi;
			Eps1inF += Eps0in*hFi;

			//Eps1out_tem = Alf1_new*hTem + Teta*(MC.Alf1M-MC.Alf1A)*hFi;
			//Eps1in_tem = Alf1_new*hTem + Teta*(MC.Alf1M-MC.Alf1A)*hFi;


			Sig1+=hSig1;
			S1+=hS1;
			Sig2+=hSig2;
			S2+=hS2;
			Sig1out = Sig1+S1;
			Sig1in = Sig1-S1;
			Sig2in = Sig2+S2;
			Sig2out = Sig2-S2;

			if (Tem >= MC.Af && (Tem - MC.Af)<=hTem) {
				TWSM = TWSM - EpsElastic;//Eps1out;////only for Alf1=0 and Alf2=0!
			}

			Tem = Tem_new;
			D1_old = D1_new;
			Alf1_old = Alf1_new;
			//EpsElastic = (Eps1out - Eps2out)/2.0;
			EpsElastic = ((MC.l_n/MC.rad)*(MC.h1+MC.h2-MC.rad*Eps2out)-MC.r_l+Eps1out*((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l))/((MC.l_n/MC.rad)*(2*MC.rad+MC.h1+MC.h2+Eps2out*MC.rad)-MC.r_l+Eps1out*((MC.l_n/MC.rad)*(MC.rad+MC.h1+MC.h2)-MC.r_l)) ;
			fprintf(res, "\n%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg", 
							Tem, M, Sig1*1e-6, S1*1e-6, Sig2*1e-6,S2*1e-6,Sig1out*1e-6,Sig1in*1e-6,Sig2in*1e-6,Sig2out*1e-6,
							Eps1out*1e2,Eps1in*1e2,Eps2in*1e2,Eps2out*1e2,Fi,6.0*M/MC.h/MC.h*1e-6,EpsElastic*1e2, MC.h2*(Eps1out-Eps1in)+MC.h1*(Eps2out-Eps2in), Eps1out_tem*1e2, Eps1in_tem*1e2);
		}//while (Tem < MC.Tmax)

			fprintf(res_sm, "\t%lg\t%lg\n", TWSM*1e2, KrBimet);
			fclose(res);
			fclose(res_sm);
			printf("\n\nExperiment finished.");
		//	pause();

			return 0;
		}
