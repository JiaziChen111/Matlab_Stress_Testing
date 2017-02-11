function [sigmult,hmat,penalt] = mhcov(p,vs)

% evaluate hessian at posterior mode

para = zeros(length(p),1);

for i = 1:length(p)
    para(i) = p(i).estval;
end

[hmat,errorij] = hessn('hessn_fcn',para,vs);

X = -hmat;

[u,s,v] = svd(X,0);

rankHHM = sum(vs.pmaskinv);

invHHMdet = 1;

for i = 1:length(p)
    if i > rankHHM
        s(i,i) = 0;
    else
        s(i,i) = 1.0/s(i,i);
        invHHMdet = invHHMdet*s(i,i);
    end
end

invHHM = u*s*u';
sigmult = u*sqrt(s);
penalt = (rankHHM/2)*log(2*pi) + 0.5*log(invHHMdet);

disp('Determinant of variance matrix');
disp(invHHMdet);
disp('s.e. at posterior mode');
disp(sqrt(diag(invHHM)));
disp('Post Mode Penalty');
disp(penalt);









