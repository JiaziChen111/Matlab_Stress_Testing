% Forecast
%
% Keith Sill
% August 2008
%

clear;

% Model setup file is run first
[M,YY,histm,lobs,nti,steps] = FcstSetup();

% Get subsample draws from DSGE model estimation
dname   = 'c:/output/matlab/26-Jan-2011/smpl-1984Q1-2010Q3-block-';
nblocks = 100;
nsim    = 10000;
nskip   = 100;
ndraws  = (nsim*nblocks)/nskip;
npara   = 40;
pmat    = SubsampleDsgeDraws(dname,nblocks,nsim,nskip,npara);

pmat    = pmat(0.2*ndraws:end,:);
mpara   = mean(pmat);
mpara(19) = 2.0;

[TT,QQ,RR,HH,DD,ZZ,VV,RC] = sysmat(mpara);


%% Since we forecastper-capita real GDP growth, we unwind it to real GDP 
% growth by adjusting for population growth --- ndot.
%steps = 12;

% population growth adjustments. We can take these from the MA simulation
ndot = 0.25*ones(steps,1);
ndot(1) = 0.25;
ndot(2) = 0.25;
ndot(3) = 0.25;
ndot(4) = 0.25;
ndot(5) = 0.25;
ndot(6) = 0.25;
ndot(7) = 0.25;
ndot(8) = 0.25;
ndot(9) = 0.25;
ndot(10) = 0.25;

% Initialize storage matricies. **1 is for shock uncertainty only case. **2
% is for shock + parameter uncertainty
[nr,nc] = size(pmat);

yfcst2    = zeros(nr,steps);
cfcst2    = zeros(nr,steps);
ifcst2    = zeros(nr,steps);
hfcst2    = zeros(nr,steps);
pfcst2    = zeros(nr,steps);
Rfcst2    = zeros(nr,steps);


% Simulate

tic;

% Generate state observations from Kalman filter
% [smat,xvar,ypred,evec,lambda,retcode] = genrstates(mpara,YY,M);
% state = smat(end,:);

for i = 1 : nr
    
	shocks = randn(steps,7);
    
    % simulations for shock uncertainty
	%fmat1 = pdrawsDSGE(mpara,state,shocks,steps);
    % Hold funds rate at zero counterfactual
    %fmat1 = pdrawsDSGEalt(mpara,state,shocks,steps);
    
    % simulations for shock and parameter uncertainty
    para_draw  = pmat(i,:);
    
    para_draw(19) = 2.0;
    
    [smat,xvar,ypred,evec,lambda,retcode] = genrstates(mpara,YY,M);
    state_draw = smat(end,:);
    
    fmat2 = pdrawsDSGE(para_draw,state_draw,shocks,steps);
    % Hold funds rate at zero counterfactual
    %fmat2 = pdrawsDSGEalt(para_draw,state_draw,shocks,steps);

    yfcst2(i,:) = fmat2(1,:);
	cfcst2(i,:) = fmat2(2,:);
	ifcst2(i,:) = fmat2(3,:);
	hfcst2(i,:) = fmat2(4,:);
	pfcst2(i,:) = fmat2(6,:);
	Rfcst2(i,:) = fmat2(7,:);

end
%%
%
% Mean forecast and coverage probabilities for the two simulation cases
%
% Real GDP growth

myfcst2 = mean(yfcst2);
pyfcst2 = hpdint(yfcst2,0.68);
ymat2   = 4 * [pyfcst2(1,:)'+ndot, myfcst2'+ndot, pyfcst2(2,:)'+ndot];
ymat2   = 400*(exp(ymat2/400)-1);

% Real Consumption growth

mcfcst2 = mean(cfcst2);
pcfcst2 = hpdint(cfcst2,0.68);
cmat2   = 4 * [pcfcst2(1,:)'+ndot, mcfcst2'+ndot, pcfcst2(2,:)'+ndot];
cmat2   = 400*(exp(cmat2/400)-1);

% Real Investment growth

mifcst2 = mean(ifcst2);
pifcst2 = hpdint(ifcst2,0.68);
imat2   = 4 * [pifcst2(1,:)'+ndot, mifcst2'+ndot, pifcst2(2,:)'+ndot];
imat2   = 400*(exp(imat2/400)-1);

% Aggregate Hours

mhfcst2 = mean(hfcst2);
phfcst2 = hpdint(hfcst2,0.68);
hmat2   = [phfcst2(1,:)', mhfcst2', phfcst2(2,:)'];


% Core PCE Inflation

mpfcst2 = mean(pfcst2);
ppfcst2 = hpdint(pfcst2,0.68);
pmat2   = 4 * ([ppfcst2(1,:)', mpfcst2', ppfcst2(2,:)'] );
pmat2   = 400*(exp(pmat2/400)-1);

% Fed funds rate

mRfcst2 = mean(Rfcst2);
pRfcst2 = hpdint(Rfcst2,0.68);
Rmat2   = ([pRfcst2(1,:)', mRfcst2', pRfcst2(2,:)'] );

toc;

nti = seqa(lobs,0.25,steps);
xlswrite('yfcst.xls',[nti,ymat2])
xlswrite('cfcst.xls',[nti,cmat2])
xlswrite('ifcst.xls',[nti,imat2])
xlswrite('hfcst.xls',[nti,hmat2])
xlswrite('pfcst.xls',[nti,pmat2])
xlswrite('Rfcst.xls',[nti,Rmat2])
