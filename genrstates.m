function [smat,smatm,ypred,evec,lambda,retcode] = genrstates(para,YY,M)
    
%/** initialize **/
[nobs,ny] = size(YY);

% /**********************************************************
% **      System matrices for State-space
% **********************************************************/
[TT,QQ,RR,HH,DD,ZZ,VV,RC]=sysmat(para);

[ns,neps] = size(RR);
Rmat = (RR'*RR)\RR';

% /**********************************************************
% ** Check determinacy 
% **********************************************************/
if (RC(1) == 1) && (RC(2)==1);
   % determinacy 
   retcode(1) = 0;

elseif (RC(1) == 1) && (RC(2)==0) 
   % indeterminacy 
   retcode(1) = 1;
   return;

else
   % no equilibrium exists, numerical problems 
   retcode(1) = RC(1);
   return;

end

% /**********************************************************
% **      Kalman filter generates state matrix
% **********************************************************/
[nstate,nx] = size(TT);
smat = zeros(nobs,ns);
ypred = zeros(nobs,ny);
evec  = zeros(nobs,neps);

% /** initial mean and variance for the state vector **/
At = zeros(nstate,1);
Pt = dlyap(TT,RR*QQ*RR');
Pt0 = Pt;

%/** filter loop **/

for t = 1:nobs

   At1 = At;
   Pt1 = Pt;

   %/** Forecasting **/
   alphahat = TT*At1;
   Phat     = TT*Pt1*TT' + RR*QQ*RR';
   Phat     = .5*(Phat+Phat');

   yhat     = ZZ*alphahat + DD;
   nut      = YY(t,:) - yhat';
   Ft       = ZZ*Phat*ZZ' + HH + ZZ*RR*VV + (ZZ*RR*VV)';
   Ft       = .5*(Ft+Ft');
   
   %/** Updating **/
   Fti = inv(Ft);
   At = alphahat + Phat*ZZ'*Fti*nut';
   Pt = Phat - (Phat*ZZ'+RR*VV)*Fti*(Phat*ZZ'+RR*VV)';

   smat(t,:) = At';
   ypred(t,:) = (DD + ZZ * At)';
   
   if t > 1
       %Reps = (smat(t,:)' - TT*smat(t-1,:)');
       %evec(t,:) = (Rmat*Reps)';
       evec(t,:) = RR\(smat(t,:)'-TT*smat(t-1,:)');
       
   end
   
   
   
end

smatm = smat * M;
lambda = diag(Pt0)'*M;



