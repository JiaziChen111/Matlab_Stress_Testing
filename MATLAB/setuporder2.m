function [KD,KDD]=setuporder2(F,xss,x,xl,e)
% F: vector of symbolic expressions defining equation system
% xss: numeric vector giving steady state value of x
% x: symbolic vector of names of current variables in system
% xl: symbolic vector of names of lagged variables in system
% e: symbolic vector of names of error terms in equations
neq=size(F,1);
nx=size(x,1);
ne=size(e,1);
nv=nx*2+ne;
scf1=jacobian(F,[x;xl;e]);
ncf1=double(subs(scf1,[x; xl; e],[xss;xss;zeros(ne,1)]));
scf2=cell(neq,1);
ncf2=cell(neq,1);
for iq=1:neq
    scf2(iq)={jacobian(scf1(iq,:),[x;xl;e])};
    ncf2(iq)=double(subs(scf2(iq),[x; xl; e],[xss;xss;zeros(ne,1)]));
end
ncf2=reshape(double([ncf2{:}]),nv,nv,neq);
ncf2=permute(ncf2,[3 1 2]);
ix0=1:nx;
ix1=ix0+nx;
ie=2*nx+(1:ne);
KD{1}=ncf1(:,ix0);
KD{2}=ncf1(:,ix1);
KD{3}=ncf1(:,ie);
KDD{1,1}=ncf2(:,ix0,ix0);
KDD{1,2}=ncf2(:,ix0,ix1);
KDD{1,3}=ncf2(:,ix0,ie);
KDD{2,2}=ncf2(:,ix1,ix1):
KDD{2,3}=ncf2(:,ix1,ie);
KDD{3,3}=ncf2(:,ie,ie);
