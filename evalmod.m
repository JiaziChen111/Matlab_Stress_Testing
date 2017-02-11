% /********************************************************************
% ** Likelihood Evaluation using Kalman Filter
% **
% ** {loglh,retcode,obserror,obsvar} = evalmod(para,YY,msel);
% ** {loglh,retcode,obserror,obsvar} = evalmod3(para,YY,msel);
% **
% ** #include ****_mod.src,  to run this procedure
% **
% ** Last updated: 06/16/2005 
% ********************************************************************/
% 

function [rloglh,retcode,obserror,obsvar] = evalmod(para,YY)

% /** likelihood conditional on first tcond observations **/
% /** to make it comparable to VAR **/
tcond = 0;      

%/** initialize **/
[nobs,ny] = size(YY);
% nobs = rows(YY);
% ny   = cols(YY);

loglhzero = -1E8;

loglh     = 0;
retcode   = 0;
obserror  = zeros(nobs,ny);
obsvar    = zeros(nobs,ny);


% /**********************************************************
% **      System matrices for State-space
% **********************************************************/
[TT,QQ,RR,HH,DD,ZZ,VV,RC]=sysmat(para);



% /**********************************************************
% ** Check determinacy 
% **********************************************************/
if (RC(1) == 1) && (RC(2)==1);
   %/* determinacy */
   retcode(1) = 0;

elseif (RC(1) == 1) && (RC(2)==0) 
   %/* indeterminacy */
   retcode(1) = 1;
   rloglh = loglhzero;
   return;

else
   %/* no equilibrium exists, numerical problems */
   retcode(1) = RC(1);
   rloglh = loglhzero;
   return;

end

% /**********************************************************
% **      Likelihood evaluation by Kalman filter
% **********************************************************/
[nstate,nx] = size(TT);

% /** initial mean and variance for the state vector **/
At = zeros(nstate,1);
Pt = dlyap(TT,RR*QQ*RR');

%/** filter loop **/
t = 1;
while t <= nobs

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
   
   %/** to make likelihood comparable to VAR(tcond), condition on first tcond observations **/
   if t > tcond;
    loglh    =  loglh -0.5*ny*log(2*pi) - 0.5*log(det(Ft)) - 0.5*nut*inv(Ft)*nut';
   end

   %/** Updating **/
   At = alphahat + Phat*ZZ'*inv(Ft)*nut';
   Pt = Phat - (Phat*ZZ'+RR*VV)*inv(Ft)*(Phat*ZZ'+RR*VV)';

   obserror(t,:) = nut;
   obsvar(t,:)   = diag(Ft)';

   t = t+1;
end

if isnan(real(loglh))==1;
   loglh = loglhzero;
end

rloglh = real(loglh);







