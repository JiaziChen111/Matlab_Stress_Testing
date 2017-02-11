function f = minout(p,g,v)

%load optim_out.mat
for i = 1:length(p)
    paraest(i) = p(i).estval;
end

[lnpy,retcode,obserror,obsvar] = evalmod(paraest,v.data);
lnprio  = priodens( paraest, v.pmean, v.pstdd, v.pshape);

for i = 1:length(p)   
    gradmod(i) = g(i)/ numgrad('invtrans',paraest(i),v.trspec);
end

disp('');
disp('===========================================');
disp(sprintf('Posterior Mode     %g',lnpy+lnprio));
disp(sprintf('Likelihood at Mode %g',lnpy));
disp(sprintf('Prior at Mode      %g',lnprio));
disp('');
disp(sprintf('%15s%15s%15s%15s%15s%15s','Parameter','Estimate','Start-V','Prior Mean','Prior Std','Gradient'));
for i = 1:length(p)
    disp(sprintf('%15s%15.5g%15.5g%15.5g%15.5g%15.5g',p(i).name,p(i).estval,p(i).sval,v.pmean(i),v.pstdd(i),gradmod(i)));
end

% 
% /********************************************************
% * Compute Residuals
% ********************************************************/
% 
% /* moments */
% 
    resid1 = obserror;
    [nr,nc] = size(obserror);
    resid1 = resid1(2:nr,:);
    [nr1,nc1] = size(resid1);
    
    resid1_m1 = mean(resid1);
    resid1_m2 = mean( (resid1 - kron(resid1_m1,ones(nr1,1))).^2 )';
    resid1_m3 = mean( (resid1 - kron(resid1_m1,ones(nr1,1))).^3 )';
    resid1_m4 = mean( (resid1 - kron(resid1_m1,ones(nr1,1))).^4 )';

    resid2 = obserror./sqrt(obsvar);
    [nr2,nc2] = size(resid2);
    resid2 = resid2(2:nr2,:);
    [nr2,nc2] = size(resid2);
    resid2_m1 = mean(resid2);
    resid2_m2 = mean( (resid2 - kron(resid2_m1,ones(nr2,1))).^2 )';
    resid2_m3 = mean( (resid2 - kron(resid2_m1,ones(nr2,1))).^3 )';
    resid2_m4 = mean( (resid2 - kron(resid2_m1,ones(nr2,1))).^4 )';


    nresid = nr1;
% 
% 
% /* autocorrelation */
% 
    ny = nc;
    resid1_ar = zeros(1,ny);
    resid2_ar = zeros(1,ny);

    for ind = 1:ny
       resid1_ar(1,ind) = inv(resid1(1:nresid-1,ind)'*resid1(1:nresid-1,ind))...
                                *resid1(1:nresid-1,ind)'*resid1(2:nresid,ind);
       resid2_ar(1,ind) = inv(resid2(1:nresid-1,ind)'*resid2(1:nresid-1,ind))...
                                *resid2(1:nresid-1,ind)'*resid2(2:nresid,ind);
    end
% 
% 
% /** print on screen **/
% if _verbose;
    disp('Raw Residuals');
    disp('-------------');
    disp(sprintf('Mean: %g\t', resid1_m1));
    disp(sprintf('Stdd: %g\t', sqrt(resid1_m2)));
    disp(sprintf('Skew: %g\t', resid1_m3./resid1_m2.^(3/2)));
    disp(sprintf('Kurt: %g\t', resid1_m4./(3*resid1_m2.^2)));
    disp(sprintf('Corr: %g\t', resid1_ar));
    disp('');
    disp('Standardized Residuals');
    disp('----------------------');
    disp(sprintf('Mean: %g\t', resid2_m1));
    disp(sprintf('Stdd: %g\t', sqrt(resid2_m2)));
    disp(sprintf('Skew: %g\t', resid2_m3./resid2_m2.^(3/2)));
    disp(sprintf('Kurt: %g\t', resid2_m4./(3*resid2_m2.^2)));
    disp(sprintf('Corr: %g\t', resid2_ar));
 