% ******************************************************************************
%  PRTR.G
%  sets parameter transformation type and arguments
% 
%  upadted:  05/31/2005
% ******************************************************************************/
function [para] = prtr()
% 
%    Each row has the following specification:
% 
%     tr~a~b~c
% 
%     tr parameter transformation type
%         0: no transformation needed
%         1: [a,b] -> [-1,1] -> [-inf,inf] by (1/c)*c*z/sqrt(1-c*z^2)
%         2: [0,inf] -> [-inf,inf] by b + (1/c)*ln(para[i]-a);
%     a  transformation argument a (usually lower bound)
%     b  transformation argument b (usually upper bound)
%     c  transformation argument c
% 
       
  alp           = [1,	1E-5,	.99,	1]; 	
  zeta_p        = [1,	1E-5,	.99,	1];
  iota_p        = [1,	0,	1.0,	1];	
  del           = [1,	1E-5,	.99,	1]; 
  ups           = [2,	0,	10,	1];
  Bigphi        = [2,	0,	5,	1];	
  s2            = [2,	0,	20,	1];	
  h             = [1,	1E-5,	.99,	1]; 	
  a2            = [2,	0,	0.99,   1];	
  nu_l          = [2,	1E-5,	10,	1];	
  nu_m          = [2,	1E-5,	100,    1];	
  zeta_w        = [1,	1E-5,	.99,	1];	
  iota_w        = [1,	0,	1,	1];	
  law           = [2,	1E-5,	50,	1];	
  rstar         = [2,	1E-5,	10,	1];	
  psi1          = [2,	1E-5,	10,	1];	
  psi2          = [2,	1E-5,	10,	1];	
  rho_r         = [1,	1E-5,	.99,	1];	
  pistar        = [0,	-1,	10,	0];	
    
  gam           = [2,	1E-6,	10,	1];	
  wadj          = [0,	0,	10,	0];	
  chi           = [2,	1E-6,	10,	1];	
  laf           = [2,	1E-5,	50,	1];	
  gstar         = [1,	1E-5,	.99,	1];	
  Ladj          = [0,	1E-6,	5000,   0];	
    
  rho_z         = [1,	0,	.99,	1];	
  rho_phi       = [1,	1E-5,	.99,	1];	
  rho_chi       = [1,	1E-5,	.99,	1];	
  rho_laf       = [1,	1E-5,	.99,	1];	
  rho_mu        = [1,	1E-5,	.99,	1];	
  rho_b         = [1,	1E-5,	.99,	1];	
  rho_g         = [1,	1E-5,	.99,	1];	
    
  sig_z         = [2,	1E-7,	100,	1];	
  sig_phi       = [2,	1E-7,	100,	1];	
  sig_chi       = [2,	1E-7,	100,	1];	
  sig_laf       = [2,	1E-7,	10000,	1];	
  sig_mu        = [2,	1E-7,	100,	1];	
  sig_b         = [2,	1E-7,	100,	1];	
  sig_g         = [2,	1E-7,	100,	1];	
  sig_r         = [2,	1E-7,	100,	1];	
    
  para = [alp ; zeta_p ; iota_p ; del ; ups ; Bigphi ; s2 ; h ; a2 ; nu_l ; nu_m ;...
      zeta_w ; iota_w ; law ;rstar ; psi1 ; psi2 ; rho_r ; pistar ; gam ; wadj ;...
      chi ; laf ; gstar ; Ladj ; rho_z ; rho_phi ;rho_chi ; rho_laf ; rho_mu ;...
      rho_b ; rho_g ; sig_z ; sig_phi ; sig_chi ; sig_laf ; sig_mu ; sig_b ;...
      sig_g ; sig_r ];
