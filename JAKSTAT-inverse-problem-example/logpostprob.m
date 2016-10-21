function [lpp PropSoln] ...
    = logpostprob(dataArray,ddefn,ddehist,N,priorpars,parsprop,tp,PropSoln)

%% Obtain the likelihood vector for each state

if  nargin == 8
    
    [logprobs PropSoln] ...
        = lik_dde(dataArray,ddefn,ddehist,N,parsprop,PropSoln);
    
elseif nargin == 7
    
    [logprobs PropSoln] ...
        = lik_dde(dataArray,ddefn,ddehist,N,parsprop);
    
end

%% Compute the log posterior probability

try
    
    lpp = sum(tp*cellfun(@sum,logprobs))+ priorlpp(priorpars,parsprop);
    
    if ~isreal(lpp)
        
        lpp         = -Inf;
        PropSoln    = NaN;
        
    end
    
catch
    
    lpp             = -Inf;
    PropSoln        = NaN;
    
end

end