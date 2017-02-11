% Set up data to estimate the DSGE model.
% Keith Sill
% May 2009
%
% The data is assumed to have been downloaded as a tab-delimited file from
% the FRED database in St Louis. The datalist in Fred has to be ordered as
% it is here (ie, the columns have to be right)
%
% Sample size is that of the imported Fred data set.
%
clear;
% Read in Quarterly data
[qvars,qdate,qraw] = xlsread('DSGE_Estim.xls','Quarterly');

% Read in Monthly data
[mvars,mdata,mraw] = xlsread('DSGE_Estim.xls','Monthly');

% Convert monthly data to quarterly by taking monthly average

[nobs,nvars] = size(mvars);

for i = 3:3:nobs
    
    j = floor(i/3);
    mqvar(j,:) = (mvars(i,:)+mvars(i-1,:)+mvars(i-2,:))/3;
end

% Data manipulations
%
coe     = qvars(:,1);
gdp     = qvars(:,2);
gdpctpi = qvars(:,3);
gpdi    = qvars(:,4);
hoanbs  = qvars(:,5);
jcxfe   = qvars(:,6);
pcdg    = qvars(:,7);
pcec    = qvars(:,8);
pcectpi = qvars(:,9);
pcesv   = qvars(:,10);
pcnd    = qvars(:,11);
 
cnp16ov = mqvar(:,1);
ff      = mqvar(:,2);
unrate  = mqvar(:,3);

[r,c] = size(gdp);

% levels transformations
ry = gdp*10^9 ./ (cnp16ov*10^3) ./ (gdpctpi/100);
rc = (pcec-pcdg)*10^9 ./ (cnp16ov*10^3) ./ (gdpctpi/100);
ri = (pcdg+gpdi)*10^9 ./ (cnp16ov*10^3) ./ (gdpctpi/100);
ah = hoanbs./cnp16ov * 4.8608*10^5;
nw = coe*10^9./(cnp16ov*10^3.*ah); 

% growth rates and logs. First obs of dataset is now 1954Q4
dry = 100 * log(ry(2:end)./ry(1:end-1));
drc = 100 * log(rc(2:end)./rc(1:end-1));
dri = 100 * log(ri(2:end)./ri(1:end-1));
hrs = log(ah(2:end));
dnw = 100 * log(nw(2:end)./nw(1:end-1));
inf = 100 * log(gdpctpi(2:end)./gdpctpi(1:end-1));
ff  = ff(2:end);

% Auxiliary variable transformations
unrate = unrate(2:end);
pcinf  = 100*log(pcectpi(2:end)./pcectpi(1:end-1));
cpcinf = 100*log(jcxfe(2:end)./jcxfe(1:end-1));

ti = seqa(1984.0,0.25,r-1);
qdata = [ ti, dry, drc, dri, hrs, dnw, inf, ff, pcinf, cpcinf, unrate ];

% data for plots and shock decomposition graphics
dpop = 400*log(cnp16ov(2:end)./cnp16ov(1:end-1));
rgdp = gdp./gdpctpi;
drgdp = 400*log(rgdp(2:end)./rgdp(1:end-1));
ninf  = 4*inf;
plotdata = [ti,dpop,drgdp,ninf,ff];

% write to a csv file
if exist('usModelData.csv','file');
    delete('usModelData.csv');
end
xlswrite('usModelData.xls',qdata);
xlswrite('graphData.xls',plotdata);

% summary statistics
disp(sprintf('var       mean         std       min       max'))
disp(sprintf('date  %10.3f  %10.3f  %10.3f  %10.3f',mean(ti),std(ti),min(ti),max(ti)))
disp(sprintf('y     %10.3f  %10.3f  %10.3f  %10.3f',mean(dry),std(dry),min(dry),max(dry)))
disp(sprintf('c     %10.3f  %10.3f  %10.3f  %10.3f',mean(drc),std(drc),min(drc),max(drc)))
disp(sprintf('i     %10.3f  %10.3f  %10.3f  %10.3f',mean(dri),std(dri),min(dri),max(dri)))
disp(sprintf('h     %10.3f  %10.3f  %10.3f  %10.3f',mean(hrs),std(hrs),min(hrs),max(hrs)))
disp(sprintf('w     %10.3f  %10.3f  %10.3f  %10.3f',mean(dnw),std(dnw),min(dnw),max(dnw)))
disp(sprintf('pi    %10.3f  %10.3f  %10.3f  %10.3f',mean(inf),std(inf),min(inf),max(inf)))
disp(sprintf('ff    %10.3f  %10.3f  %10.3f  %10.3f',mean(ff),std(ff),min(ff),max(ff)))
disp(sprintf('pcinf %10.3f  %10.3f  %10.3f  %10.3f',mean(pcinf),std(pcinf),min(pcinf),max(pcinf)))
disp(sprintf('cpcinf%10.3f  %10.3f  %10.3f  %10.3f',mean(cpcinf),std(cpcinf),min(cpcinf),max(cpcinf)))
disp(sprintf('unrate%10.3f  %10.3f  %10.3f  %10.3f',mean(unrate),std(unrate),min(unrate),max(unrate)))

figure;
subplot(2,2,1);
plot(ti',dry);
title('dy');

subplot(2,2,2);
plot(ti',drc);
title('dc');

subplot(2,2,3);
plot(ti',dri);
title('di');

subplot(2,2,4);
plot(ti',hrs);
title('hrs');

figure;
subplot(2,2,1);
plot(ti',dnw);
title('dw');

subplot(2,2,2);
plot(ti',inf);
title('inf');

subplot(2,2,3);
plot(ti',ff);
title('ff');


