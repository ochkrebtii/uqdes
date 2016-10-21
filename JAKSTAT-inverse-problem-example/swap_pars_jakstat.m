function [newpars,newprobsoln,new_log_r_bot,swap_counter] ...
    = swap_pars_jakstat(old2pars,old2propsoln,old2logbot,...
                        old2tp,swap_counter,prop_pars,PropSoln,...
                        dataArray,priorpars,ddefn,ddehist,N)

% This will check if the parameters should be swapped between the two
% chains.  Note that although it would make sense to keep updating
% DEfd_chain, it can't be done in the parfor loop.  I think that this is
% because I need to save it into cell array locations as it is currently coded.
% The parfor loop makes a big difference so it's worth keeping in this form.
% This function works if both chains use smooth approximations to the ODE
% or if one or both of them uses the actual ODE solution.
% Use parallel_temper_tp if using parallel tempering in which case this
% is the divisor of the log top of the acceptance ratio.

% old2pars        = {chainpars{chains_swap(1)}(:,iter),chainpars{chains_swap(2)}(:,iter)};
% old2propsoln    = {PropSoln{chains_swap(1)},PropSoln{chains_swap(2)}};
% old2logbot      = [log_r_bot{chains_swap(1)}(iter),log_r_bot{chains_swap(2)}(iter)];
% old2tp          = tp(chains_swap);
% swap_counter    = 0;


try
    
    for tt=1:2
        
        [log_r_top(tt) PropSoln_temp{tt}] ...
            = logpostprob(dataArray,ddefn,ddehist,N,priorpars,...
            prop_pars{tt},old2tp(tt),PropSoln{tt});
        
    end
    
    if isinf(old2logbot(2)) && ~isinf(old2logbot(1))
        
        % if the chain with more flexibility (cooler temp) is okay but the higher
        % temp chain is stuck with isinf(log_r_bot)==1, then we should assume that
        % the lower temp chain is obligated to pass on its parameter value to the
        % higher temp chain.
        
        newpars                         = prop_pars;
        newprobsoln                     = PropSoln;
        swap_counter                    = swap_counter + 1;
        new_log_r_bot                   = {log_r_top(1),log_r_top(2)};
        
    else
        
        if exp(sum(log_r_top)-sum(old2logbot))>rand
            
            newpars                     = prop_pars;
            newprobsoln                 = PropSoln;
            swap_counter                = swap_counter + 1;
            new_log_r_bot               = {log_r_top(1),log_r_top(2)};
            
        else
            
            newpars                     = old2pars;
            newprobsoln                 = old2propsoln;
            new_log_r_bot               = {old2logbot(1),old2logbot(2)};
            
        end
        
    end
    
    new_log_r_bot                       = cell2mat(new_log_r_bot);
    
catch
    
    newpars                             = old2pars;
    newprobsoln                         = old2propsoln;
    new_log_r_bot                       = old2logbot;
    
    warning('problem in swap')
    
end


end

