function [TT,QQ,RR,HH,DD,ZZ,VV,RC] = sysmat(para)

  % Set names to parameters
  
  alp	 = para(1);
  del	 = para(4);	
  ups	 = exp(para(5)/100);
  Bigphi = para(6);
   h	 = para(8);
   nu_m	 = para(11);
   bet 	 = exp(-para(15)/400);
  pistar = exp(para(19)/400);
  gam	 = para(20)/400;
  chi	 = para(22);
  laf	 = para(23);
  gstar	 = 1+para(24);
  Ladj	 = para(25);
  sigz	 = para(33);
  sigphi = para(34);
  siglaf = para(36);
  sigmu	 = para(37);
  sigb	 = para(38);
  sigg	 = para(39);
  sigR	 = para(40);
  
  % Implicit parameter definitions
  Lstar	= 1;
  rstar = (1/bet)*exp(gam)*ups^(alp/(1-alp)); 
  rkstar= (1/bet)*exp(gam)*ups^(1/(1-alp))-(1-del);
  omegastar= (alp^(alp)*(1-alp)^(1-alp)*rkstar^(-alp)/(1+laf))^(1/(1-alp));
  kstar = (alp/(1-alp))*omegastar*Lstar/rkstar;
  kbarstar= kstar*exp(gam)*ups^(1/(1-alp));
  istokbarst= 1-((1-del)/(exp(gam)*ups^(1/(1-alp))));
  istar = kbarstar*istokbarst;
  ystar = (kstar^alp)*(Lstar^(1-alp))-Bigphi;
  
  cstar = (ystar/gstar)-istar;
  xistar = (1/cstar)*((1/(1-h*exp(-gam)*ups^(-alp/(1-alp))))-(h*bet/(exp(gam)*ups^(alp/(1-alp))-h)));
  Rstarn = pistar*rstar;
  wstar = ( 1/(1+laf) * (alp^alp) * (1-alp)^(1-alp) * rkstar^(-alp) )^(1/(1-alp));
  mstar = ( chi * Rstarn/(1-Rstarn) * xistar^(-1) )^(1 / nu_m);

  if ystar <= 0
    disp('ystar is negative!');
    disp(alp);
    disp(bet);
    disp(kstar);
    disp(Lstar);
  end
  
  Rstarn = pistar*rstar;
  
  %
  %      Solve DSGE model using GENSYS
  %

  [T1,TC,T0,RC,GAM0,GAM1] = dsgesolv(para);

  %
  %      System matrices
  %      Initialization
  %

  % observation indicies

  eq_y  = 1;
  eq_c  = 2;
  eq_i  = 3;
  eq_h  = 4;
  eq_w  = 5;
  eq_pi = 6;
  %eq_m  = 7;
  eq_R  = 7;
  %eq_cy = 8;
  %eq_iy = 9;
  %eq_wy = 10;
  %eq_my = 12;

  % number of observation variables
  ny = 7;

  % variable indicies (match with dsgesolv indicies for T1)
  v_y 	= 1;
  v_c 	= 2;
  v_i 	= 3;
  v_L 	= 6;
  v_m 	= 7;
  v_pi 	= 9;
  v_R 	= 10;
  v_w 	= 13;
  v_z   = 17;

  % Shock indicies (neps)
  e_z	= 1;
  e_phi	= 2;
  e_mu	= 3;
  e_b	= 4;
  e_g	= 5;
  e_laf	= 6;
  e_R   = 7;

  %
  %       Transition equations with augmentation
  %
  %       z(t) = T1*z(t-1) + T0*e(t)
  %       x(t) = [z(t)',...)'
  %
  %       x(t) = TT*x(t-1) + RR*e(t)
  %       e(t) ~ iid N(0,QQ)
  %

  [rz,nz] = size(T1);
  %nz  = cols(T1);
  [rep,nep] = size(T0);
  %nep = cols(T0);

  % augmentation: lagged variables 
  v_yl = nz+1;
  v_cl = nz+2;
  v_il = nz+3;
  v_wl = nz+4;
  v_ml = nz+5;

  TT = zeros(nz+5,nz+5);

  TT(1:nz,1:nz)  = T1;
  TT(v_yl,v_y)   = 1;
  TT(v_cl,v_c)   = 1;
  TT(v_il,v_i)   = 1;
  TT(v_wl,v_w)   = 1;
  TT(v_ml,v_m)   = 1;

  QQ = zeros(nep,nep);
  QQ(e_z,e_z)     = sigz^2;
  QQ(e_phi,e_phi) = sigphi^2;
  QQ(e_mu,e_mu)   = sigmu^2;
  QQ(e_b,e_b)     = sigb^2;
  QQ(e_g,e_g)     = sigg^2;
  QQ(e_laf,e_laf) = siglaf^2;
  QQ(e_R,e_R)     = sigR^2;

  RR = [T0; zeros(5,nep)];

  %
  %  Measurement equations
  %
  %  y(t) = DD + ZZ*x(t) + u(t)
  %  u(t) ~ N(0,HH)
  %  
  %  cov(e(t),u(t)) = VV
  %
  
  [rs,nstate] = size(TT);

  DD = zeros(ny,1);           % demeaned 
  
  DD(eq_y,1)  = 100*(gam+(alp*log(ups)/(1-alp)));
  DD(eq_c,1)  = 100*(gam+(alp*log(ups)/(1-alp)));
  DD(eq_i,1)  = 100*(gam+(alp*log(ups)/(1-alp)));
  DD(eq_h,1)  = log(Ladj) + log(Lstar);
  DD(eq_w,1)  = 100*(gam + log(pistar) + (alp*log(ups)/(1-alp)));
  DD(eq_pi,1) = 100*log(pistar);
  %DD(eq_m,1)  = (log(pistar)+gam+ alp*log(ups)/(1-alp));
  DD(eq_R,1)  = 400*log(Rstarn);
  
  %DD(eq_cy,1) = 100*log(cstar/ystar);
  %DD(eq_iy,1) = 100*log(istar/ystar);
  %DD(eq_wy,1) = 100*log(wstar/ystar);
  %DD(eq_my,1) = 100*log(mstar/ystar);
  
  ZZ = zeros(ny,nstate);
 
  
  ZZ(eq_y,v_y)  =  1;
  ZZ(eq_y,v_yl) = -1;
  ZZ(eq_y,v_z)  = 1;
  
  ZZ(eq_c,v_c)  = 1;
  ZZ(eq_c,v_cl) = -1;
  ZZ(eq_c,v_z)  = 1;
  
  ZZ(eq_i,v_i)  = 1;
  ZZ(eq_i,v_il) = -1;
  ZZ(eq_i,v_z)  = 1;
  
  ZZ(eq_h,v_L) = 1/100;
  
  ZZ(eq_w,v_w)  = 1;
  ZZ(eq_w,v_wl) = -1;
  ZZ(eq_w,v_z)  = 1;
  ZZ(eq_w,v_pi) = 1;
  
  ZZ(eq_pi,v_pi) = 1;
  
  %ZZ(eq_m,v_m)  = 1;
  %ZZ(eq_m,v_ml) = -1;
  %ZZ(eq_m,v_z)  = 1;
  
  ZZ(eq_R,v_R) = 4;
  
  %ZZ(eq_cy,v_c) = 1;
  %ZZ(eq_cy,v_y) = -1;
  
%   ZZ(eq_iy,v_i) = 1;
%   ZZ(eq_iy,v_y) = -1;
%   
%   ZZ(eq_wy,v_w) = 1;
%   ZZ(eq_wy,v_y) = -1;
  
  %ZZ(eq_my,v_m) = 1;
  %ZZ(eq_my,v_y) = -1;
  
  HH = zeros(ny,ny);  % w/o measurement error 
  
  VV = zeros(nep,ny);
  
end
  

