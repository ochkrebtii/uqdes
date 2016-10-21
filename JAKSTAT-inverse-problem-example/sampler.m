function [old_pars old_probsoln log_r_bot accepts] ...
    = sampler(old_pars,old_probsoln,log_r_bot,stepvar,tpn,accepts,dataArray,ddefn,ddehist,N,priorpars)

% this function draws a new sample, accepts/rejects it, and outputs
% the new parameters and the log of the denominator of the rejection
% ratio

%% Define sampling blocks of parameters


% old_pars          = chainpars{tc}(:,iter);
% old_probsoln      = PropSoln{tc};
% log_r_bot         = log_r_bot{tc}(iter);
% stepvar           = stepvar_chain{tc};
% tpn               = tp(tc);
% accepts           = accepts{tc}(k,:);


sample.blocind{1}         = 1:8;                          % index of block parameters
sample.bloctarget{1}      = [0.18,0.28];                  % target interval for acceptance rate
sample.resolve{1}         = 1;
sample.posprop{1}         = 0;

% sample.blocind{2}         = 5:6;                          % index of block parameters
% sample.bloctarget{2}      = [0.23,0.33];                  % target interval for acceptance rate
% sample.resolve{2}         = 0;

% sample.blocind{2}         = 7;                            % index of block parameterst
% sample.bloctarget{2}      = [0.33,0.43];                  % target interval for acceptance rate
% sample.resolve{2}         = 1;

% sample.blocind{3}         = 8;                            % index of block parameterst
% sample.bloctarget{3}      = [0.33,0.43];                  % target interval for acceptance rate
% sample.resolve{3}         = 1;

sample.blocind{2}         = 9:10;                         % index of block parameterst
sample.bloctarget{2}      = [0.34,0.44];                  % target interval for acceptance rate
sample.resolve{2}         = 1;
sample.posprop{2}         = 1;

% sample.blocind{6}         = 11:12;                         % index of block parameterst
% sample.bloctarget{6}      = [0.23,0.33];                  % target interval for acceptance rate
% sample.resolve{6}         = 1;

sample.blocnum            = length(sample.blocind);  % instead of parameters, return number of sampled blocks



if (nargin == 0) % if calling sample() then only return information about sampling blocks
    
    old_pars                    = sample;
    old_probsoln                = [];
    log_r_bot                   = [];
    accepts                     = [];
    return
    
end

for ind = 1:sample.blocnum
    
    new_pars                                = old_pars;
    
    if sample.posprop{ind} == 0
        new_pars(sample.blocind{ind})           = mvnrnd(old_pars(sample.blocind{ind}),stepvar(sample.blocind{ind},sample.blocind{ind}))';
    else
        new_pars(sample.blocind{ind})           = abs(mvnrnd(old_pars(sample.blocind{ind}),stepvar(sample.blocind{ind},sample.blocind{ind}))');
    end
    
    [log_r_top new_probsoln] ...
        = logpostprob(dataArray,ddefn,ddehist,N,priorpars,new_pars,tpn);
    
    % calculate rejection ratio
    if sample.posprop{ind} == 0
        logrho  = log_r_top - log_r_bot;
    else
        logrho  = log_r_top + mvncdf(old_pars(sample.blocind{ind}),zeros(size(old_pars(sample.blocind{ind}))),stepvar(sample.blocind{ind},sample.blocind{ind})) ...
            - log_r_bot - mvncdf(new_pars(sample.blocind{ind}),zeros(size(old_pars(sample.blocind{ind}))),stepvar(sample.blocind{ind},sample.blocind{ind}));
    end
    
    % accept/reject
    if unifrnd(0,1) < exp(logrho)                                           % accept
        
        old_pars                            = new_pars;                     % update new parameters
        accepts(ind)                        = accepts(ind) + 1;             % update acceptance count
        log_r_bot                           = log_r_top;                    % new value of log denominator to output
        old_probsoln                        = new_probsoln;
        
    end
    
end


end