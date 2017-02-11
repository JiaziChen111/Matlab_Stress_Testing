function [retcode,Atmat,Ptmat] = kfilt(para,YY)

% Run the Kalman Filter, return the sequences {S(t|t)} and {P(t|t)}

[nobs,ny] = size(YY);

retcode = 0;

% Get state space representation 
[TT,QQ,RR,HH,DD,ZZ,VV,RC] = sysmat(para);
%TT = xlsread('../gaussv4/ttmat.xls');
%RR = xlsread('../gaussv4/rrmat.xls');

% /**********************************************************
% ** Check determinacy 
% **********************************************************/
if (RC(1) == 1) && (RC(2)==1);
   %/* determinacy */
   retcode(1) = 0;

elseif (RC(1) == 1) && (RC(2)==0) 
   %/* indeterminacy */
   retcode(1) = 1;
   return;

else
   %/* no equilibrium exists, numerical problems */
   retcode(1) = RC(1);
   return;

end

[nstate,nx] = size(TT);

% /** initial mean and variance for the state vector **/
At = zeros(nstate,1);
Pt = dlyap(TT,RR*QQ*RR');

Atmat = zeros(nobs,nstate);
Ptmat = zeros(nobs,nstate*nstate);

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

   %/** Updating **/
   At = alphahat + Phat*ZZ'/ Ft * nut';
   Pt = Phat - (Phat*ZZ'+RR*VV) / Ft * (Phat*ZZ'+RR*VV)';
   
   Atmat(t,:) = At';
   Ptmat(t,:) = reshape(Pt,1,nstate*nstate);

   t = t+1;
end