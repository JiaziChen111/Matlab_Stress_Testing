function mparamat = SubsampleDsgeDraws(dname,nblocks,nsim,nskip,npara)
%
% Subsample draws from the Metropolis Hastings run on the DSGE model
%

paramat = zeros(nblocks*nsim,npara);

for i = 1:nblocks

    clear fname parasim;

    fname = strcat(dname,num2str(i));
    load(fname,'parasim','rej');
    ridx1 = (i-1)*nsim + 1;
    ridx2 = i*nsim;
    paramat(ridx1:ridx2,:) = parasim;

end

ndraws   = (nblocks*nsim)/nskip;
mparamat = zeros(ndraws,npara);


for i = 1:ndraws
    mparamat(i,:) = paramat(i*nskip,:);
end

