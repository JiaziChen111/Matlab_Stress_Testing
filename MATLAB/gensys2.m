function [FD,FDD,M11,M2,C,q,zs,zu,v,gev,eu]=gensys2(KD,KDD,c,Pi,omega,pick,div)
%function [FD,FDD,M11,M2,C,q,zs,zu,v,gev,eu]=gensys2(KD,KDD,c,Pi,omega,pick,div)
%------- Including Kollman repair to constant term --------------------------
% System originally in the form K(z(t),z(t-1),eps(t),eta(t))=0, t=0,...,\infty, 
% where Et[eta(t+1)]=Et[eps(t+1)]=0 for t\ge 0. Its second order expansion is, for equation j,
%  KD{1}(j,:)dz(t)= KD{2}(j,:)dz(t-1) + KD{3,:}(j)eps(t) + Pi_j eta(t) ...
%     + .5*(dz(t)'*KDD{1,1}(j,:,:)*dz(t) + 2*dz(t)'*KDD{1,2}(j,:,:)*dz(t-1)) ...
%     + 2*dz(t)'*KDD{1,3}(j,:,:)*eps(t) + dz(t-1)'*KDD{2,2}(j,:,:)*dz(t-1) ...
%     +2*dz(t-1)'*KDD{2,3}(j,:,:)*eps(t)+eps(t)'*KDD{3,3}(j,:,:)*eps(t)) + Pi*eta +C
% (note that the above notation is not quite right as matlab code.  KDD{2,3}{j,:,:) is not
% an ordinary 2x2 matrix until it has been "squeezed") 
% It is assumed that the eta terms enter linearly and thus generate no 
% second-order terms. 
% sigma^2*omega=var(eps(t))
% pick*dz(t) is the proposed state vector.  If pick is present, then if possible, the solution 
% will be constructed with dy(t)=pick*dz(t).  Otherwise a default method is used to find a
% state vector.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
%  existence only with not-s.c. eps; eu=[-2,-2] for coincident zeros. 
%  eu(1) incremented by 10 for existence problems with 2nd order solution.
% Solution is y(t)=FD{1}*z(t-1)+FD{2}*eps(t)+FD{3}*sigma^2 ...
%                +.5*prodt(prodt(FDD{1,1},w(t-1),2,1),w(t-1),2,1)
%                 +.5*prodt(prodt(FDD{1,2},w(t-1),2,1),eps(t),2,1)
%                 +.5*permute(prodt(FDD(1,2},w(t-1),2,1),eps(t),2,1),[1,3,2]
%                 +.5*prodt(FDD(2,2),eps(t),2,1),eps(t),2,1)
%                  +C(six)
%             x(t)=.5*prodt(prodt(M11,y(t),2,1),y(t),2,1)+M2*sigma^2+C(usix)
% The FDD{j,k} terms for k<j are not calculated, as they are just permutations of the 
% corresponding k,j terms.  To get back to the original z vector, use [y(t);x(t)]=[zs;zu]z(t)
% or [pick;zu]z(t), according to whether pick has been supplied and used or not.
% Use gensys.m for first-order solutions that go backward and forward, allowing
% serial correlation in eps.  This program requires serially uncorrelated eps, but
% otherwise provides all the analysis in gensys.m, plus second-order analysis.
% KD (3x1), KDD(3x3), FD(3x1), and FDD(2x2) are cell arrays whose elements are numerical
% arrays.  The elements of KDD and FDD are three-subscript arrays.
% By Christopher A. Sims
% 7/27/00
eu=[0;0];
realsmall=1e-8;
fixdiv=(nargin==7);
n=size(KD{1},1);
[a b q z v]=qz(KD{1},KD{2});
if ~fixdiv, div=1e10; end
nunstab=0;
zxz=0;
for i=1:n
%------------------div calc------------
   if ~fixdiv
      if abs(a(i,i)) > 0
         divhat=abs(b(i,i))/abs(a(i,i));
         if 1+realsmall<divhat & divhat<div
            div=.5*(1+divhat);
         end
      end
   end
%----------------------------------------
   nunstab=nunstab+(abs(b(i,i))>div*abs(a(i,i)));
   if abs(a(i,i))<realsmall & abs(b(i,i))<realsmall
      zxz=1;
   end
end
%div 
%nunstab
if ~zxz
   [a b q z,v]=qzdiv(div,a,b,q,z,v);
end
gev=[diag(a) diag(b)];
if zxz
   %disp('Coincident zeros.  Indeterminacy and/or nonexistence.')
   eu=[-2;-2];
   return
end
six=1:n-nunstab;
usix=n-nunstab+1:n;
q1=q(six,:);
q2=q(usix,:);
zs=z(:,six)';
zu=z(:,usix)';
a2=a(usix,usix);
b2=b(usix,usix);
etawt=q2*Pi;
zwt=q2*KD{3};
[ueta,deta,veta]=svd(etawt);
md=min(size(deta));
bigev=find(diag(deta(1:md,1:md))>realsmall);
ueta=ueta(:,bigev);
veta=veta(:,bigev);
deta=deta(bigev,bigev);
[uz,dz,vz]=svd(zwt);
md=min(size(dz));
bigev=find(diag(dz(1:md,1:md))>realsmall);
uz=uz(:,bigev);
vz=vz(:,bigev);
dz=dz(bigev,bigev);
if isempty(bigev)
   exist=1;
else
   exist=norm(uz-ueta*ueta'*uz) < realsmall*n;
end
if ~isempty(bigev)
   zwtx0=b2\zwt;
   zwtx=zwtx0;
   M=b2\a2;
   for i=2:nunstab
      zwtx=[M*zwtx zwtx0];
   end
   zwtx=b2*zwtx;
   [ux,dx,vx]=svd(zwtx);
   md=min(size(dx));
   bigev=find(diag(dx(1:md,1:md))>realsmall);
   ux=ux(:,bigev);
   vx=vx(:,bigev);
   dx=dx(bigev,bigev);
   existx=norm(ux-ueta*ueta'*ux) < realsmall*n;
else
   existx=1;
end
%----------------------------------------------------
% Note that existence and uniqueness are not just matters of comparing
% numbers of roots and numbers of endogenous errors.  These counts are
% reported below because usually they point to the source of the problem.
%------------------------------------------------------
[ueta1,deta1,veta1]=svd(q1*Pi);
md=min(size(deta1));
bigev=find(diag(deta1(1:md,1:md))>realsmall);
ueta1=ueta1(:,bigev);
veta1=veta1(:,bigev);
deta1=deta1(bigev,bigev);
if existx | nunstab==0
   %disp('solution exists');
   eu(1)=1;
else
    if exist
        %disp('solution exists for unforecastable z only');
        eu(1)=-1;
    %else
        %fprintf(1,'No solution.  %d unstable roots. %d endog errors.\n',nunstab,size(ueta1,2));
    end
    %disp('Generalized eigenvalues')
   %disp(gev);
   %md=abs(diag(a))>realsmall;
   %ev=diag(md.*diag(a)+(1-md).*diag(b))\ev;
   %disp(ev)
%   return;
end
if isempty(veta1)
   unique=1;
else
   unique=norm(veta1-veta*veta'*veta1)<realsmall*n;
end
if unique
   %disp('solution unique');
   eu(2)=1;
else
   fprintf(1,'Indeterminacy.  %d loose endog errors.\n',size(veta1,2)-size(veta,2));
   %disp('Generalized eigenvalues')
   %disp(gev);
   %md=abs(diag(a))>realsmall;
   %ev=diag(md.*diag(a)+(1-md).*diag(b))\ev;
   %disp(ev)
%   return;
end
tmat = [a(six,six)\[eye(n-nunstab) -(ueta*(deta\veta')*veta1*deta1*ueta1')'];...
                  zeros(nunstab,n-nunstab) eye(nunstab)];
%--------- go through order2guts without worrying about pick ------------------
%if nargin>5 % so pick was specified
%   %if size(pick)==size(zs) & cond(pick*zs')>realsmall
%   % line above makes no sense.  cond() always returns >1.  
%   if size(pick)==size(zs)& norm(pick*zs')>realsmall
%%      tmat=[pick*z;zeros(nunstab,n-nunstab) eye(nunstab)]*tmat; 
%%      z=inv([pick*z;zu]);
%      tmat(six,:)=(pick*zs')*tmat(six,:);
%      z=inv([pick;zu]);
%   else
%      disp('pick matrix does not select a true state.  Default state used.')
%   end
%end
%----------
BD=KD;
BDD=KDD;
for j=1:3
   BD{j}=tmat*q*KD{j};
   for k=j:3
      BDD{j,k}(:)=tmat*q*KDD{j,k}(:,:);
   end
end
for j=1:2
   BD{j}=BD{j}*z;
   BDD{j,3}=permute(prodt(BDD{j,3},z,2,1,3),[1 3 2]);
   for k=j:2
      BDD{j,k}=prodt(prodt(BDD{j,k},z,2,1,3),z,2,1,3);
   end
end
%-----------------------
C=tmat*q*c;
C(usix)=(BD{1}(usix,usix)-BD{2}(usix,usix))\C(usix);

% CORRECTION RK 26/3/2001
%C(six)=C(six)+(BD{2}(six,usix)-BD{1}(six,usix))*C(usix);
% The correction above was needed only because the previous version eliminated lagged 
% x spuriously.  
F1=BD{2}(six,:);
F2=BD{3}(six,:);
if max(abs(gev(six,1)).\abs(gev(six,2)))*max(abs(gev(usix,2)).\abs(gev(usix,1))) >= 1.0
   eu(1)=eu(1)+10;
   rc=1;
   M11=[];M2=[];FDD{1,1}=[];FDD{1,2}=[];FDD{2,2}=[];F3=[];
else
   [M11,M2,FDD,F3,rc]=order2guts(BD,BDD,F1,F2,omega) ;
end
FD{1}=F1;FD{2}=F2;FD{3}=F3;
