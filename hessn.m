% Compute Hessian

function [hessian,errorij] = hessn(fcn,xc,v)

n = length(xc);
ndx = 6;
for i = 1:ndx
    dx(i) = 5 + 2*(i-1);
end
dx = exp(-dx);

dxscale = ones(n,1);

hessian = zeros(n);
hessdiag = zeros(ndx,1);

fx = feval(fcn,xc,v);

fmt = '%5.0f ';

disp(' ');
disp(' ');
disp('Computing Hessian..');
disp('--------------------------------------------------');
disp(sprintf(['  diagonal elements     :   ' fmt ],n));

% Diagonal elements first
for i = 1:n
    
    for j = 1:ndx
        paradx = xc;
        parady = xc;
        paradx(i) = paradx(i) + dx(j)*dxscale(i);
        parady(i) = parady(i) - dx(j)*dxscale(i);
        fdx = feval(fcn,paradx,v);
        fdy = feval(fcn,parady,v);
        hessdiag(j) = -(2*fx - fdx -fdy)/(dx(j)*dxscale(i))^2;
    end
    
    hessian(i,i) = -0.5*(hessdiag(3)+hessdiag(4));
    disp(sprintf(['                            ' fmt ],i));
end

% Off-diagonal elements next
disp(' ');
disp(sprintf(['  off-diagonal elements :   ' fmt ],n*(n-1)/2));

errorij = [];
k = 1;

for seli = 1:n
    for selj = seli+1:n
        
        for i = 1:ndx
            paradx = xc;
            parady = xc;
            paradx(seli) = paradx(seli) + dx(i)*dxscale(seli);
            parady(selj) = parady(selj) - dx(i)*dxscale(selj);
            paradxdy = paradx;
            paradxdy(selj) = paradxdy(selj) - dx(i)*dxscale(selj);
            fdx = feval(fcn,paradx,v);
            fdy = feval(fcn,parady,v);
            fdxdy = feval(fcn,paradxdy,v);
            hessdiag(i) = -(fx - fdx - fdy + fdxdy)/(dx(i)*dx(i)*dxscale(seli)*dxscale(selj));
        end
        
        hessian(seli,selj) = -0.5*(hessdiag(3)+hessdiag(4));
        
        if hessian(seli,seli) == 0 || hessian(selj,selj) == 0
            corrij = 0;
        else
            corrij = hessian(seli,selj)/sqrt(hessian(seli,seli)*hessian(selj,selj));
        end
        
        if abs(corrij) > 0.98
            hessian(seli,selj)= 0.9*sqrt(hessian(seli,seli)*hessian(selj,selj))
            errorij = [errorij,seli,selj,corrij];
        elseif abs(corrij) < 0.005
            hessian(seli,selj) = 0.0;
        end
        
        hessian(selj,seli) = hessian(seli,selj);
        
        if mod(k,5) == 0
            disp(sprintf(['                            ' fmt ],k));
        end
        k = k+1;
    end
end

hessian = real(hessian);
