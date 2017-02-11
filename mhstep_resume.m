function mhstep_resume(sigmult,hessian,nblock,nsim,cc0,cc,para,vst,runname,blockno)

npara = max(size(para));    

% best to display some model identification stuff here?

sigscale = sigmult;
sigpropinv = hessian;

sigpropdim = sum(vst.pmaskinv);

[u,s,v] = svd(sigpropinv);
sigproplndet = 0.0;

for i = 1:npara
    if i <= sigpropdim
        s(i,i) = 1.0/s(i,i);
        sigproplndet = sigproplndet + log(s(i,i));
    end
end

paranew = para';

tic;
for n = 1:nblock
    
    parasim     = zeros(nsim,npara);
    likesim     = zeros(nsim,1);
    postsim     = zeros(nsim,1);
    rej         = zeros(nsim,1);
    propsim     = zeros(nsim,1);
    proppostsim = zeros(nsim,1);
    
    likesim(1,1) = likenew;
    postsim(1,1) = postnew;
    propsim(1,1) = propdens;
    proppostsim(1,1) = postnew;
    parasim(1,:) = paranew';
    postold = postnew;
    likeold = likenew;
    paraold = paranew;
    
    newj = 2;
    
    for j = newj:nsim
        
        paranew = paraold + cc*(sigscale*randn(npara,1));
     
        [postnew,likenew] = mhstep_fcn(paranew,vst);
    
        propdens = -0.5*sigpropdim*log(2*pi) - 0.5*sigproplndet - ...
            0.5*sigpropdim*log(cc*cc) - 0.5*(paranew-paraold)'*sigpropinv*(paranew-paraold)/(cc*cc);
    
        propsim(j,1) = propdens;
        proppostsim(j,1) = postnew;
        
        r = min(min(1,exp(postnew-postold)));
        
        if rand < r
            postsim(j,1) = postnew;
            likesim(j,1) = likenew;
            parasim(j,:) = paranew';
            
            paraold = paranew;
            postold = postnew;
            likeold = likenew;
            
        else
            likesim(j,1) = likeold;
            postsim(j,1) = postold;
            parasim(j,:) = paraold';
            rej(j) = 1;
            
        end
        
        if mod(j,nsim/10) == 0;
            disp(sprintf('block %d of %d : sim % d of %d',n,nblock,j,nsim));
        end
        
    end
    
    disp(sprintf('Rejection percent: %f\n',mean(rej)));
    disp(sprintf('Likelihood:        %f\n',mean(likesim)));
    disp(sprintf('Posterior:         %f\n',mean(postsim)));
    
    fname = strcat(runname,'-block-',num2str(n+blockno));
    
    save(fname, 'parasim', 'postsim', 'likesim', 'proppostsim', 'rej');
    
    
    
end

disp(sprintf('elapsed time is %f minutes',(toc/60)));
