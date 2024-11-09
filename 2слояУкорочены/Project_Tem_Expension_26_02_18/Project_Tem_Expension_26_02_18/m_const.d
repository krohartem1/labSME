
// Material constants for steel and TiNi elemens
  
 SY2 = 2000 //1200yield limit for steel in MPa
  SY1 = 50 //280 // initial phase yield limit for TiNi in MPa
  SY1A = 600 //600 // 800 initial yield limit for TiNi in austenitic state in MPa
      E2 = 20.0//26//23.5 //23.5 //25 //20.8 //18.0 //21.7 // Young modulus of steel in GPa in MPa
  E1A = 29.9//32.3//39//32.3 //28.3 //31.6 //14.4 //34.0// Young_A in GPa
  E1M_L= 17.5//18.6//17.6 //12.5// 25.0 //Young_M at loading in GPa
  E1M_U =17.5//18.6 //17.6 //12.5	 // Young_M at unloading in GPa
  H2 = 9  // strengthening coefficient for steel in GPa
  H1 = 4.0 //15.0 // strengthening coefficient for TiNi in GPa
  H1A = 5 //10 // strengthening coefficient for TiNi in austenitic state in GPa
  h1 = 0.01 // height  TiNi in mm
  h2 = 0.03 // height steel in mm
  h = 0.04
  b1 = 2 // width of steel element in mm
  b2 = 2 // width of TiNi element in mm
  k_b = 0.1 // ratio b1/b2
  k_h = 0.1//ratio h1/(h1+h2)
  Kr = 1.0 // recovery ratio
  r_l = 0.8//0.7 //0.95 //3.7 //3.7 //0.7 // length difference in mm
  l_n = 124.0 // length of plate in mm
  rad = 20 // radius in mm
  EpsT = 1.0 //preliminary strain in pct
  Alf1A =15.4 //0.0 // Heat expansion of TiNi, 10^-6
  Alf1M =9.2 //0.0 // Heat expansion of TiNi, 10^-6
  Alf2 = 10.0 //0.0 // Heat expansion of TiNi, 10^-6

Ms = 312//313 //K
Mf = 302//311 // K
As = 306//314 //K
Af = 316 //K
Tmin = 273//305//K
Tmax = 353//320//320 //K

Lam1 = 0.0
ETP = 15.2 // in GPa

M_t = 100 //total moment in N*m
 