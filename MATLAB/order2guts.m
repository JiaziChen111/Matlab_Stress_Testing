function [M11,M2,FDD,F3,rc]=order2guts(BD,BDD,F1,F2,omega)     
%  function [M11,M2,FDD,F3,rc]=order2guts(BD,BDD,F1,F2,omega)     
%  System originally in the form Ka(z(t),z(t-1),eps(t),eta(t))=0, t=0,...,\infty, 
%  where Et[eta(t+1)]=0 for t\ge 0. Its second order expansion is K1 dz(t)= K2 
%  dz(t-1) + K3\eps(t) + K4\eta(t) + .5*(K11*dz(t)\otimes dz(t) + K12*dz(t)\otimes 
%  dz(t-1) ... It is assumed that the \eta terms enter linearly and thus generate no 
%  second-order terms. gensys transformations convert this form  to an equation 
%  where the K's are replaced by B's and z is are replaced by w=[y,x], with the 
%  following special characteristics:  the xy block of BD{1} is zero; the yy block of 
%  BD{1} is the identity; the y row of BD{4} is zero (and hence BD{4} is unused); the 
%  xx block of BD{2} is 
%  non-singular, and BD{2}{usix,usix}\BD{1}{usix,usix} (usix indexes the xx block of BD{2}) has 
%  all its eigenvalues 
%  less than (1/div)<1 in absolute value; the yy block of BD{2} has all its eigenvalues 
%  less than div>1 in absolute value. F1 and F2 are the coefficients in the first order
%  solution y(t)=F1*w(t-1)+F2*eps(t).  omega is the covariance matrix of eps.
%  Note that BD and BDD are cell arrays, whose elements are matrices.  BDD{j,k} for k<j
%  are not used, so should be left null to save space.
%  The "state" may be interpreted either as y(t) or as s(t)=F1*w(t).  Often these have the
%  same row space, but they need not always.  When they don't, z(t) is a function of s(t-1) for
%  all t, while it cannot be written as a function of y(t-1) alone at the initial date.  
% 
[nstate,n1]=size(F1); %if nstate~=nstate2, error('F1 matrix not square'), end;
[n,n2]=size(BD{1}); if n~=n2, error('BD{1} matrix not square'),end;   
uindex=nstate+1:n; 
sindex=1:nstate;
Fs=F1(sindex,sindex);
nu=n-nstate;
nerr=size(F2,2);
L=BD{2}(uindex,uindex)\BD{1}(uindex,uindex); 
M11 = -prodt(BDD{1,2}(uindex,sindex,sindex),Fs,2,1,3); 
M11=M11+permute(M11,[1 3 2])-prodt(prodt(BDD{1,1}(uindex,sindex,sindex),Fs,2,1,3),Fs,2,1,3)...
   -BDD{2,2}(uindex,sindex,sindex);
M11(:)=BD{2}(uindex,uindex)\M11(:,:);
%------------------Doubling algorithm-------------------------
R=Fs; 
eest=1; 
Minc=zeros(size(M11)); 
crit= 10*sqrt((n-nstate)*nstate*nstate)*eps; 
dblinc=0;
while eest>crit & dblinc<200
   Minc(:)=L*M11(:,:);
   M11=M11+prodt(prodt(Minc,R,2,1),R,2,1);
   L=L*L;
   R=R*R;
   eest=sqrt(sum(sum(sum(Minc.*Minc))));
   dblinc=dblinc+1;
end
%fprintf(1,'dblinc=%d, eest=%g\n',dblinc,eest);
if dblinc>200
   fprintf(1,'Inaccuracy %g+%gi in order2guts.m\n ', [real(eest) imag(eest)])
   rc=1;
else
   rc=0;
end
worka=prodt(prodt(M11,F1,2,1,3),F1,2,1,3);
worka=reshape(-BD{1}(sindex,uindex)*worka(:,:),[nstate,n,n]);
% in non-ic version, workb=prodt(BDD{1,2}(sindex,sindex,sindex!),F1,2,1,3);
workb=prodt(BDD{1,2}(sindex,sindex,:),F1,2,1,3);
workb=workb+permute(workb,[1 3 2]);
%FDD{1,1}=prodt(prodt(BDD{1,1}(sindex,sindex,sindex),F1,2,1),F1,2,1)+worka+workb...
%   +BDD{2,2}(sindex,sindex,sindex);
% CORRECTION, APRIL 12, 2001, by Kollman 
%workc=prodt(BD{2}(sindex,uindex),M11,2,1,2);
% added in ic:---
%workc=cat(2,workc,zeros(nstate,nu,nstate));
%workc=cat(3,workc,zeros(nstate,n,nu));
%---
% Kollman correction above was compensating for the original error of not
% allowing for the possibility that y depends on lagged x as well s lagged y.
%
FDD{1,1}=prodt(prodt(BDD{1,1}(sindex,sindex,sindex),F1,2,1,3),F1,2,1,3)+worka+workb...
%  +workc + BDD{2,2}(sindex,:,:);
+ BDD{2,2}(sindex,:,:);
%----------------
worka=prodt(prodt(M11,F2,2,1,3),F2,2,1,3);
worka(:)=BD{1}(uindex,uindex)*worka(:,:);
workb=prodt(prodt(BDD{1,1}(uindex,sindex,sindex),F2,2,1,3),F2,2,1,3);
workc=prodt(BDD{1,3}(uindex,sindex,:),F2,2,1,3);
workc=workc+permute(workc,[1 3 2]);
worka=worka-workb-workc-BDD{3,3}(uindex,:,:);
M2=.5*((BD{2}(uindex,uindex)-BD{1}(uindex,uindex))\(worka(:,:)*omega(:)));
%----------------
worka=prodt(prodt(M11,F1,2,1,3),F2,2,1,3);
%pre-ic:  worka=reshape(-BD{1}(sindex,uindex)*worka(:,:),[nstate,nstate,nerr]);
worka=reshape(-BD{1}(sindex,uindex)*worka(:,:),[nstate,n,nerr]);
% line below involves a change by Kollman, plus an additional change by CAS to make
% it work when the state vector is one dimensional.
%workb=cat(2,prodt(BDD{1,2}(sindex,sindex,sindex),F2,2,1,3)...
%   +permute(prodt(BDD{1,3}(sindex,sindex,:),F1,2,1,3),[1 3 2])+BDD{2,3}(sindex,sindex,:);
workb=cat(2,prodt(BDD{1,2}(sindex,sindex,sindex),F2,2,1,3),zeros(nstate,nu,nerr))...
   +permute(prodt(BDD{1,3}(sindex,sindex,:),F1,2,1,3),[1 3 2])+BDD{2,3}(sindex,:,:);
FDD{1,2}=prodt(prodt(BDD{1,1}(sindex,sindex,sindex),F1,2,1,3),F2,2,1,3)+worka+workb;
%--------------
worka=prodt(prodt(M11,F2,2,1,3),F2,2,1,3);
worka=reshape(-BD{1}(sindex,uindex)*worka(:,:),[nstate,nerr,nerr]);
worka=worka+prodt(prodt(BDD{1,1}(sindex,sindex,sindex),F2,2,1,3),F2,2,1,3);
workb=prodt(BDD{1,3}(sindex,sindex,:),F2,2,1,3);
workb=workb+permute(workb,[1 3 2]);
FDD{2,2}=worka+workb+BDD{3,3}(sindex,:,:);
%---------------------
% ORIGINAL CODE
F3=-BD{1}(sindex,uindex)*M2;

% CORRECTION, APRIL 12, 2001, by Kollman
% F3=(-BD{1}(sindex,uindex)+BD{2}(sindex,uindex))*M2; 
