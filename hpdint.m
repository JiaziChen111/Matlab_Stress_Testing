function hpdband = hpdint(draws,percent)

[ndraws,drawdim] = size(draws);
hpdband = zeros(2,drawdim);
nwidth = floor(percent*ndraws);

for i = 1:drawdim
    
    drawcoli = draws(:,i);
    drawcoli = sort(drawcoli,'descend');
    bup = 1;
    minwidth = drawcoli(1) - drawcoli(nwidth);
    
    for j = 2:(ndraws-nwidth+1)
        newwidth = drawcoli(j) - drawcoli(j+nwidth-1);
        if newwidth < minwidth
            bup = j;
            minwidth = newwidth;
        end
    end
    
    hpdband(2,i) = drawcoli(bup);
    hpdband(1,i) = drawcoli(bup+nwidth-1);
    
end
