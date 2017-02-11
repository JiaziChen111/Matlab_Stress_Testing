function [M,YY,histm,lobs,nti,steps] = FcstSetup()

% Get estimation sample period data and set the state selection matrix M
%

% Form the reduced-state selection matrix
%
% The state variable ordering is:
% y,c,i,k,,kbar,L,m,mc,pi,R,rk,u,w,wtil,xi,xik,z,phi,mu,b,g,laf,Ec,Epi,
% Erk,Ew,Ewtil,Exi,Exik,Ei,yl,cl,il,wl,ml
%
% The 'true' state variables are c,i,kbar,R,w,z,phi,mu,b,g,laf which are
% indexed as 2,3,5,10,13,17,18,19,20,21,22
%
% note that pi is not a state when the estimated model shuts down inflation
% indexation, as we have done here.
%
M = zeros(35,11);
M(2,1)   = 1;
M(3,2)   = 1;
M(5,3)   = 1;
M(10,4)  = 1;
M(13,5)  = 1;
M(17,6)  = 1;
M(18,7)  = 1;
M(19,8)  = 1;
M(20,9)  = 1;
M(21,10) = 1;
M(22,11) = 1;

% load data to start the forecast simulations
%YY   = xlsread('usModelData.xlsx','sheet3','b119..h224');
YY = dlmread('usModelData.txt','\t',118,1);


steps = 10;
sobs = 2007.0;
lobs = 2011.25;
hobs = (lobs-sobs)/0.25+1;

nti = seqa(2007,0.25,hobs+steps);

% set up history file for plots
%[rdata,txt] = xlsread('usModelData.xlsx','a157:m263');
rdata = dlmread('usModelDataRaw.txt','\t',148,1);

dpop = 100*log(rdata(2:end,8)./rdata(1:end-1,8));
dry  = 100*log((rdata(2:end,1)./rdata(2:end,7) ./ ...
    (rdata(1:end-1,1)./rdata(1:end-1,7))));
inf  = YY(:,6); %100*log((rdata(2:end,7)./rdata(1:end-1,7)));
ffed = rdata(2:end,12);

nobs = max(size(inf));
histm = [seqa(1984,0.25,nobs),dry,inf,ffed,dpop];




end