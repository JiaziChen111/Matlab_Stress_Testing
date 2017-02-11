 % nkpm_m.m
% 
% finds the posterior mode for the linearized model using the kalman filter
% to construct the likelihood
%	$Id: nkmp_pm.m,v 1.2 2008/02/13 15:05:28 keith Exp keith $	

% Initialize the model

clear;

diary;

datename = strcat('c:/output/matlab/',date,'/');
mkdir(datename);

copyfile('nkmp_pm.m',datename);
copyfile('usModelData.xlsx',datename);
copyfile('benchmark_start.csv',datename);
copyfile('benchmark_prior.csv',datename);

% Define a structure to hold parameter names, starting values, and
% estimated values

paras = struct('name',...
    {'alp','zeta_p','iota_p','del','ups','Bigphi','s2','h','a2','nu_l',...
    'nu_m','zeta_w','iota_w','law','beta','psi1','psi2','rho_r','pistar',...
    'gam','Wadj','chi','laf','gstar','Ladj','rho_z','rho_phi','rho_chi',...
    'rho_laf','rho_mu','rho_b','rho_g','sig_z','sig_phi','sig_chi',...
    'sig_laf','sig_mu','sig_b','sig_g','sig_r'},'sval',{0},'estval',{0});

% Priors --- get the prior from csv file

[prior] = csvread('benchmark_prior.csv',16,1);
pshape   = prior(:,1);
pmean    = prior(:,2);
pstdd    = prior(:,3);
pmask    = prior(:,4);
pfix     = prior(:,5);
pmaskinv = 1-pmask;
pshape   = pshape.*pmaskinv;

% Starting values for estimation in a structure that gets passed to
% csminwel

[p0vals] = csvread('benchmark_start.csv',2,1);
for i = 1:length(p0vals)
    paras(i).sval = p0vals(i);
end
npara = max(size(p0vals));

% get the parameter transformation specification. This is used to change
% the parameter vector support to help csminwell find the max

trspec = prtr();

% Load data for estimation

smpl = '-1984Q1-2011Q1';
%YY   = xlsread('usModelData.xlsx','sheet3','b119:h224');
YY = dlmread('usModelData.txt','\t',118,1);


[nobs,ny] = size(YY);

disp(sprintf('initial observations are:'));
disp(sprintf('%10.5f',YY(1,:)));
disp(sprintf('last observations are:'));
disp(sprintf('%10.5f',YY(end,:)));
disp(sprintf('program paused, hit any key'));
pause;

runname  = strcat(datename,'smpl',smpl);
estname  = strcat(datename,'estimval',smpl,'.mat');
hessname = strcat(datename,'hessn-',smpl,'.mat');
 
% Prior and Posterior at Starting Values
% Define csminwel input, 
 
% Define a structure containing variables/parameters/data that will
% get passed to the objective function minimized by csminwel
% 
 varargin = struct('pmaskinv',pmaskinv,'pfix',pfix,'pmask',pmask,'pmean',...
     pmean,'pstdd',pstdd,'pshape',pshape,'data',YY,'trspec',trspec);
 
disp('');
disp('==============================================');
lnprio = priodens(p0vals,pmean,pstdd,pshape);
disp(sprintf('Prior Density at Starting Value is: %g',lnprio)); 
lnpy = -objfcn(invtrans(p0vals,varargin.trspec),varargin);
disp(sprintf('Posterior at Starting Value is: %g',lnpy));
%disp('Press any key to continue');
%pause;
% 
% Maximize Posterior Density
% Maximize the posterior density, calculate Hessian, display results
% 
[paraest,g,retcode] = objmin(paras,varargin);
minout(paraest,g,varargin);
save(estname,'paraest','g','retcode');

% Compute Hessian at Posterior Mode
% compute square root inverse hessian

load(estname,'paraest','g','retcode');
[sigmult,hessian,penalt] = mhcov(paraest,varargin);
save(hessname,'sigmult','hessian','penalt');
  
 % Metropolis-Hastings RW Simulation
load(estname,'paraest');
load(hessname,'sigmult','hessian','penalt');

nblock = 100;
nsim   = 1000;
cc0    = 0;
cc     = 0.1; 

% Code block to restart the MHRW if it is interrupted or wanders
% off into a bad region. If a resume, set appendmode to 1 and set the
% value for the last good block

appendmode    = 0;
lastGoodBlock = 0;

if appendmode == 1
    fname = strcat(runname,'-block-',num2str(lastGoodBlock));
    load(fname,'parasim');
    paravals = parasim(nsim,1:npara);
    % Set new starting value to last good parameter draw
    for i = 1:npara
        paraest(i).estsval = paravals(i);
    end   
end

% Run the posterior simulator
mhstep(sigmult,hessian,nblock,nsim,cc0,cc,paraest,varargin,runname,...
        appendmode,lastGoodBlock);




            
