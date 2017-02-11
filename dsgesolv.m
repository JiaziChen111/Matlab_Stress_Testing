function [T1,TC,T0,RC,GAM0,GAM1] = dsgesolv(para)

  
  % Set names to parameters

  alp	 = para(1);		
  zeta_p = para(2);		
  iota_p = para(3);		
  del	 = para(4);		
  ups	 = exp(para(5)/100);		
  Bigphi = para(6);		
  s2	 = para(7);		
  h	     = para(8);		
  a2	 = para(9);		
  nu_l	 = para(10);		
  nu_m	 = para(11);		
  zeta_w = para(12);		
  iota_w = para(13);		
  law	 = para(14);		
  bet 	 = exp(-para(15)/400);		
  psi1	 = para(16);		
  psi2	 = para(17);		
  rho_r	 = para(18);		
  pistar = exp(para(19)/400);		
  gam	 = para(20)/400;		
  wadj	 = para(21);		
  chi	 = para(22);		
  laf	 = para(23);		
  gstar	 = 1+para(24);	
  Ladj	 = para(25);		
  rho_z	 = para(26);		
  rho_phi= para(27);		
  rho_chi= para(28);		
  rho_laf= para(29);		
  rho_mu = para(30);		
  rho_b	 = para(31);		
  rho_g	 = para(32);		
  sig_z	 = para(33);		
  sig_phi= para(34);		
  sig_chi= para(35);		
  sig_laf= para(36);		
  sig_mu = para(37);		
  sig_b	 = para(38);		
  sig_g	 = para(39);		
  sig_r	 = para(40);		

  % define implicit parameters

  Lstar	= 1;
  zstar = gam+(alp/(1-alp))*log(ups); 
  rstar = (1/bet)*exp(gam)*ups^(alp/(1-alp)); 
  rkstar= (1/bet)*exp(gam)*ups^(1/(1-alp))-(1-del);
  omegastar= (alp^(alp)*(1-alp)^(1-alp)*rkstar^(-alp)/(1+laf))^(1/(1-alp));
  kstar = (alp/(1-alp))*omegastar*Lstar/rkstar;
  kbarstar= kstar*exp(gam)*ups^(1/(1-alp));
  istokbarst= 1-((1-del)/(exp(gam)*ups^(1/(1-alp))));
  istar = kbarstar*istokbarst;
  ystar = (kstar^alp)*(Lstar^(1-alp))-Bigphi;

  %if ystar <= 0
  %disp([alp,  bet, kstar,Lstar])
  %dm([ystar,Lstar,kstar,Bigphi])
  %keyboard;
  %end

  cstar = (ystar/gstar)-istar;
  xistar = (1/cstar)*((1/(1-h*exp(-gam)*ups^(-alp/(1-alp))))-(h*bet/(exp(gam)*ups^(alp/(1-alp))-h)));
  phi = Lstar^(-nu_l)*omegastar*xistar/(1+law);
  Rstarn = pistar*rstar;
  wstar = ( 1/(1+laf) * (alp^alp) * (1-alp)^(1-alp) * rkstar^(-alp) )^(1/(1-alp));
  mstar = ( chi * Rstarn/(1-Rstarn) * xistar^(-1) )^(1 / nu_m);

  % Equations
  eq_margcost 	= 1;
  eq_prsett 	= 2;
  eq_capacc 	= 3;
  eq_effcap 	= 4;
  eq_euler 	    = 5;
  eq_moneydem 	= 6;
  eq_margut 	= 7;
  eq_invfoc 	= 8;
  eq_rettocap 	= 9;
  eq_utcap 	    = 10;
  eq_optwage 	= 11;
  eq_aggwage 	= 12;
  eq_caplabrat 	= 13;
  eq_resources 	= 14;
  eq_prod 	    = 15;
  eq_taylor 	= 16;
  eq_z		    = 17;
  eq_phi	    = 18;
  eq_mu		    = 19;
  eq_b		    = 20;
  eq_g		    = 21;
  eq_laf	    = 22;
  eq_Ec		    = 23;
  eq_Epi	    = 24;
  eq_Erk	    = 25;
  eq_Ew		    = 26;
  eq_Ewtil	    = 27;
  eq_Exi	    = 28;
  eq_Exik	    = 29;
  eq_Ei		    = 30;

  neq = 30;
  
  v_y  	 = 1;
  v_c 	 = 2;
  v_i 	 = 3;
  v_k 	 = 4;
  v_kbar = 5;
  v_L 	 = 6;
  v_m 	 = 7;
  v_mc 	 = 8;
  v_pi 	 = 9;
  v_R 	 = 10;
  v_rk 	 = 11;
  v_u 	 = 12;
  v_w 	 = 13;
  v_wtil = 14;
  v_xi 	 = 15;
  v_xik  = 16;

  % Exogenous variables (equations) and exogenous shocks (if shocks are iid the # of exog equations is less than the # of shocks)
  % note what we call z_t here is really z^*_t in the paper
  v_z 	= 17;
  v_phi = 18;
  v_mu 	= 19;
  v_b 	= 20;
  v_g 	= 21;
  v_laf = 22;  
  
  % Expectation terms 
  v_Ec 	 = 23;
  v_Epi  = 24;
  v_Erk  = 25;
  v_Ew 	 = 26;
  v_Ewtil= 27;
  v_Exi  = 28;
  v_Exik = 29;
  v_Ei 	 = 30;

  neq = 30;
  
  % Shock indicies (neps)
  e_z   = 1;
  e_phi = 2;
  e_mu  = 3;
  e_b   = 4;
  e_g   = 5;
  e_laf = 6;
  e_R   = 7;
  
  neps 	= 7;

  % Expectation error indicies (neta)
  n_c	 = 1;
  n_pi	 = 2;
  n_rk	 = 3;
  n_w	 = 4;
  n_wtil = 5;
  n_xi	 = 6;
  n_xik	 = 7;
  n_i	 = 8;
  
  neta = 8;
  
  GAM0 = zeros(neq,neq);
  GAM1 = zeros(neq,neq);
  C = zeros(neq,1);
  PSI = zeros(neq,neps);
  PPI = zeros(neq,neta);
  
  % marginal cost equation (1.2.56)
  GAM0(eq_margcost,v_mc) = 1;
  GAM0(eq_margcost,v_w) = -(1-alp);
  GAM0(eq_margcost,v_rk) = -alp;
  

  % price setting (1.2.59)
  GAM0(eq_prsett,v_pi) = (1+iota_p*bet)*zeta_p;
  GAM0(eq_prsett,v_Epi) = -bet*zeta_p;					  
  GAM0(eq_prsett,v_mc) = -(1-zeta_p)*(1-bet*zeta_p);
  GAM0(eq_prsett,v_laf) = -1; % transform laf_t, old version: -laf*(1-zeta_p)*(1-bet*zeta_p)/(1+laf);
  
  GAM1(eq_prsett,v_pi) = iota_p*zeta_p;
  
  
  % capital accumulation	 (1.2.65)
  GAM0(eq_capacc,v_kbar) = 1;
  GAM0(eq_capacc,v_z) = (1-istokbarst);  
  GAM0(eq_capacc,v_i) = -istokbarst;
  GAM0(eq_capacc,v_mu) = -istokbarst*(1+bet)*exp(2*zstar)*s2; % transformed mu_t; old version -istokbarst
  
  GAM1(eq_capacc,v_kbar) = (1-istokbarst);
  
  
  % effective capital  (1.2.64)
  GAM0(eq_effcap,v_k) = 1;
  GAM0(eq_effcap,v_u) = -1;
  GAM0(eq_effcap,v_z) = 1;
  
  GAM1(eq_effcap,v_kbar) = 1;


  % Euler equation (1.2.61)
  GAM0(eq_euler,v_xi) = (exp(zstar)-h*bet)*(exp(zstar)-h); 
  GAM0(eq_euler,v_b) = -(exp(2*zstar)+h^2*bet);               % old version: -(exp(zstar)-h)*exp(zstar);
  GAM0(eq_euler,v_b) = GAM0(eq_euler,v_b) + bet*h*rho_b*exp(-zstar)*(exp(2*zstar)+h^2*bet); 
  
  GAM0(eq_euler,v_z) = h*exp(zstar);
  GAM0(eq_euler,v_z) = GAM0(eq_euler,v_z)-h*bet*exp(zstar)*rho_z;   % this is Ez_t+1 
  GAM0(eq_euler,v_c) = (exp(2*zstar)+h^2*bet);
  GAM0(eq_euler,v_Ec) = -h*bet*exp(zstar);
  
  GAM1(eq_euler,v_c) = h*exp(zstar);
  
  % money demand	(1.2.62) -- no money demand shock, b shock is dropped too?
  GAM0(eq_moneydem,v_m) = nu_m;
  GAM0(eq_moneydem,v_xi)	= 1;
  GAM0(eq_moneydem,v_R) = 1/(Rstarn-1);
  
  
  % marginal utility (1.2.63)
  GAM0(eq_margut,v_xi) = 1;
  GAM0(eq_margut,v_Exi) = -1;
  GAM0(eq_margut,v_z) = rho_z;   % this is Ez_t+1 
  GAM0(eq_margut,v_R) = -1;
  GAM0(eq_margut,v_Epi) = 1;


  % investment FOC (1.2.66)
  GAM0(eq_invfoc,v_xik) = exp(-2*zstar)/s2;
  GAM0(eq_invfoc,v_mu) = (1+bet);             % new_mu is transformed as (1+bet)*exp(2*zstar)*s2*old_mu
  GAM0(eq_invfoc,v_xi) = -exp(-2*zstar)/s2;
  GAM0(eq_invfoc,v_i) = -(1+bet);
  GAM0(eq_invfoc,v_Ei) = bet;
  GAM0(eq_invfoc,v_z) = -1;
  GAM0(eq_invfoc,v_z) = GAM0(eq_invfoc,v_z) + bet*rho_z;  % this is Ez_t+1 
  
  GAM1(eq_invfoc,v_i) = -1;
  
  % return to capital (1.2.67)
  GAM0(eq_rettocap,v_xik) = 1;
  GAM0(eq_rettocap,v_Exi) = -rkstar/(rkstar+1-del);
  GAM0(eq_rettocap,v_Erk) = -rkstar/(rkstar+1-del);
  GAM0(eq_rettocap,v_Exik) = -(1-del)/(rkstar+1-del);
  GAM0(eq_rettocap,v_z) 	= rho_z; % this is Ez_t+1 
  
  
  % capital utilization (1.2.68)
  GAM0(eq_utcap,v_u) = a2;
  GAM0(eq_utcap,v_rk) = -rkstar;
  
  % optimal wage (1.2.69)
  GAM0(eq_optwage,v_wtil) = 1+nu_l*(1+law)/law;
  GAM0(eq_optwage,v_w) = 1+zeta_w*bet*nu_l*(1+law)/law;
  GAM0(eq_optwage,v_Ewtil) = -zeta_w*bet*(1+nu_l*(1+law)/law);
  GAM0(eq_optwage,v_Ew) = -zeta_w*bet*(1+nu_l*(1+law)/law);
  GAM0(eq_optwage,v_b) = -(1-zeta_w*bet)*(exp(2*zstar)+h^2*bet)*exp(-zstar)/(exp(zstar)-h);
  GAM0(eq_optwage,v_phi) = -1; % transformed phi_t; old version: -(1-zeta_w*bet);
  
  GAM0(eq_optwage,v_L) = -(1-zeta_w*bet)*nu_l;
  GAM0(eq_optwage,v_xi) = (1-zeta_w*bet);
  GAM0(eq_optwage,v_Epi) = -zeta_w*bet*(1+nu_l*(1+law)/law);
  GAM0(eq_optwage,v_z) = -zeta_w*bet*rho_z*(1+nu_l*(1+law)/law);
  GAM0(eq_optwage,v_z) = GAM0(eq_optwage,v_z) + zeta_w*bet*iota_w*(1+nu_l*(1+law)/law);
  GAM0(eq_optwage,v_pi) 	= zeta_w*bet*iota_w*(1+nu_l*(1+law)/law);
  
  % aggregate wage evolution (1.2.70)
  GAM0(eq_aggwage,v_w)    = 1;
  GAM0(eq_aggwage,v_pi)   = 1;
  GAM0(eq_aggwage,v_z)    = 1;
  GAM0(eq_aggwage,v_wtil) = -(1-zeta_w)/zeta_w;
  
  GAM1(eq_aggwage,v_w)    = 1;
  GAM1(eq_aggwage,v_pi)   = iota_w;
  GAM1(eq_aggwage,v_z) = iota_w;
  
  % capital labor ratio (1.2.60)
  GAM0(eq_caplabrat,v_k) = 1;
  GAM0(eq_caplabrat,v_L) = -1;
  GAM0(eq_caplabrat,v_w) = -1;
  GAM0(eq_caplabrat,v_rk)	= 1;
  
  % aggregate resources (1.2.71)
  GAM0(eq_resources,v_y) = 1;
  GAM0(eq_resources,v_c) = -cstar/(cstar+istar);
  GAM0(eq_resources,v_i) = -istar/(cstar+istar);
  GAM0(eq_resources,v_g) = -1;
  GAM0(eq_resources,v_u) = -rkstar*kstar/(cstar+istar);
  
  
  % production function (1.2.72)
  GAM0(eq_prod,v_y) = 1;
  GAM0(eq_prod,v_k) = -alp*(ystar+Bigphi)/ystar;
  GAM0(eq_prod,v_L) = -(1-alp)*(ystar+Bigphi)/ystar;
  
  
  % Taylor rule (1.2.73)
  GAM0(eq_taylor,v_R) = 1;
  GAM0(eq_taylor,v_pi) = -(1-rho_r)*psi1;
  GAM0(eq_taylor,v_y) = -(1-rho_r)*psi2;     % targeting output 
  GAM1(eq_taylor,v_R) = rho_r;
  
  PSI(eq_taylor,e_R) = 1;
  
  % z shock
  GAM0(eq_z,v_z) = 1;
  GAM1(eq_z,v_z) = rho_z;
  PSI(eq_z,e_z) = 1;
  
  % phi shock
  GAM0(eq_phi,v_phi) = 1;
  GAM1(eq_phi,v_phi) = rho_phi;
  PSI(eq_phi,e_phi) = 1;
  
  % mu shock
  GAM0(eq_mu,v_mu) = 1;
  GAM1(eq_mu,v_mu) = rho_mu;
  PSI(eq_mu,e_mu) = 1;
  
  % b shock
  GAM0(eq_b,v_b) = 1;
  GAM1(eq_b,v_b) = rho_b;
  PSI(eq_b,e_b) = 1;
  
  % g shock
  GAM0(eq_g,v_g) = 1;
  GAM1(eq_g,v_g) = rho_g;
  PSI(eq_g,e_g) = 1;
  
  % laf shock
  GAM0(eq_laf,v_laf) = 1;
  GAM1(eq_laf,v_laf) = rho_laf; 
  PSI(eq_laf,e_laf) = 1;
  
  %
  % expectational equations
  %
  
  % E_c
  GAM0(eq_Ec,v_c)  = 1;
  GAM1(eq_Ec,v_Ec) = 1;
  PPI(eq_Ec,n_c)   = 1;
  
  % E_pi
  GAM0(eq_Epi,v_pi)  = 1;
  GAM1(eq_Epi,v_Epi) = 1;
  PPI(eq_Epi,n_pi)   = 1;
  
  % E_rk
  GAM0(eq_Erk,v_rk)  = 1;
  GAM1(eq_Erk,v_Erk) = 1;
  PPI(eq_Erk,n_rk)   = 1;
  
  % E_w
  GAM0(eq_Ew,v_w)  = 1;
  GAM1(eq_Ew,v_Ew) = 1;
  PPI(eq_Ew,n_w)   = 1;
  
  % E_wtil
  GAM0(eq_Ewtil,v_wtil)  = 1;
  GAM1(eq_Ewtil,v_Ewtil) = 1;
  PPI(eq_Ewtil,n_wtil)   = 1;
  
  % E_xi
  GAM0(eq_Exi,v_xi)  = 1;
  GAM1(eq_Exi,v_Exi) = 1;
  PPI(eq_Exi,n_xi)   = 1;
  
  % E_xik
  GAM0(eq_Exik,v_xik)  = 1;
  GAM1(eq_Exik,v_Exik) = 1;
  PPI(eq_Exik,n_xik)   = 1;
  
  % E_i
  GAM0(eq_Ei,v_i)  = 1;
  GAM1(eq_Ei,v_Ei) = 1;
  PPI(eq_Ei,n_i)   = 1;
  
  %
  %      QZ(generalized Schur) decomposition by GENSYS
  %
  
  [T1,TC,T0,M,TZ,TY,GEV,RC] = gensys(GAM0,GAM1,C,PSI,PPI,1);
  
end
