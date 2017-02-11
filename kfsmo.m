function [sAtmat,sPtmat] = kfsmo(Atmat,Ptmat,para)

% Kalman smoother

[nobs,nstate] = size(Atmat);

[TT,QQ,RR,HH,DD,ZZ,VV,RC] = sysmat(para);
%TT = xlsread('../gaussv4/ttmat.xls');
%RR = xlsread('../gaussv4/rrmat.xls');

ahat = Atmat(nobs,:)';
phat = reshape(Ptmat(nobs,:),nstate,nstate);

sAtmat(nobs,:) = Atmat(nobs,:);
sPtmat(nobs,:) = Ptmat(nobs,:);

t = 2;

while t <= nobs
    
    at = Atmat(nobs-t+1,:)';
    pt = reshape(Ptmat(nobs-t+1,:),nstate,nstate);
    
    pnxt = TT * pt * TT' + RR * QQ * RR';
    
    [u_svd,s_svd,v_svd] = svd(pnxt);
    for i = 1:nstate
        if i > rank(pnxt);
            s_svd(i,i) = 0;
        else
            s_svd(i,i) = 1/s_svd(i,i);
        end
    end
    
    invpnxt = u_svd*s_svd*u_svd';
    
    pstr = pt * TT' * invpnxt;
    
    ahat = at + pstr * (ahat - TT * at);
    phat = pt + pstr * (phat - pnxt) * pstr';
    
    sAtmat(nobs-t+1,:) = ahat';
    sPtmat(nobs-t+1,:) = reshape(phat,1,nstate*nstate);
    
    t = t+1;
    
end