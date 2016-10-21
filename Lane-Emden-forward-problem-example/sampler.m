function [new_pars,new_X,log_r_bot,accepts] ...
    = sampler(Nknots,tspan,tgrid,old_pars,old_X,log_r_bot,prior_pars,...
                time,data,accepts,stepvar,odefn,LogPostProb,tp)

% this function draws a new sample, accepts/rejects it, and outputs
% the new parameters and the log of the denominator of the rejection
% ratio
% 

% tc = 1;
% tp=tp(tc)
% old_X = X{tc};
% old_pars = chainpars{tc}(:,iter)
% log_r_bot = log_r_bot{tc}(iter)
% prior_pars=prior_pars
% accepts=accepts{tc}(k,:)
% stepvar=stepvar_chain{tc}


%ReLpp_if_zero = 0;
% if zero, this variable indicates that none of the MH steps were accepted


%% propose new DE parameters while keeing tausqs and smoothing parameters fixed
temp_pars = [old_pars(1:2);normrnd(old_pars(3),stepvar);old_pars(4:end)];
%temp_pars = [old_pars(1:2);mvnrnd(old_pars(3),stepvar);old_pars(4:end)];

% calculate rejection ratio
[log_r_top t_X]  = LogPostProb(Nknots,tspan,tgrid,temp_pars,old_X,time,data,prior_pars,odefn,tp);
logrho          = log_r_top - log_r_bot;               % rejection ratio

% accept/reject
if (unifrnd(0,1) < exp(logrho))              % accept
    new_pars = temp_pars;               % update new parameters
    accepts(1) = accepts(1) + 1;      % update acceptance count
    log_r_bot = log_r_top;            % new value of log denominator to output
    new_X = t_X;
    %ReLpp_if_zero = 1;
else                                % reject
    new_pars = old_pars;        % keep old parameters, log_r_bot stays the same
    new_X = old_X;
end

end


