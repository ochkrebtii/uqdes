function [logprobs PropSoln] ...
    = lik_dde(dataArray,ddefn,ddehist,N,parsprop,PropSoln)

%% Sample solution realization

% need only one solution realization
nsolves                             = 1;   

% extract model and auxiliary parameters from array
%[modpars,tau,histvec,alpha,lambda,alphaepo,lambdaepo]  = extractparms(parsprop);   

[modpars,tau,histvec,alpha,lambda]  = extractparms(parsprop);   

pars.model                          = vertcat(modpars,tau);
%pars.model                          = [modpars,tau,alphaepo,lambdaepo];
%pars.solver                         = [alpha,lambda];
pars.EpoRAobs                       = dataArray.EpoRAobs;
pars.EpoRAtime                      = dataArray.EpoRAtime;

if nargin == 5 % if an old solution is not passed into the function
    
    PropSoln = uqdes_dde(dataArray.tspan,nsolves,N,tau,...
            lambda,alpha,ddefn,ddehist,pars,histvec,...
            dataArray.evalgrid);
  
    if isnan(sum(sum(PropSoln)))
    
        logprobs            = -Inf;
        PropSoln            = NaN;
        return
    
    end

    
end


%% Compute the likelihood centered at solution realization

logprobs                = cell(1,dataArray.nstates);

temp                    = obstransform(PropSoln,modpars(5:6));

for state = 1:dataArray.nstates
    
    diff                = dataArray.obs{state} - temp(dataArray.inds{state},state);
    logprobs{state}     = -(diff.^2)'*(1./dataArray.sigsq{state});

    if isnan(logprobs{state})
    
    logprobs            = -Inf;
    PropSoln            = NaN;
    
    return
    
    end
    
end



end