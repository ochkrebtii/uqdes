function [newpars,new_X,log_r_bot,swap_counter]  = ...
    swap_pars_emden(Nknots,tspan,tgrid,old_pars,old_X,log_r_bot,...
    prop_pars,prop_X,prior_pars,time,data,odefn,swap_counter,LogPostProb,tp)

% This will check if the parameters should be swapped between the two
% chains.  Note that although it would make sense to keep updating
% DEfd_chain, it can't be done in the parfor loop.  I think that this is
% because I need to save it into cell array locations as it is currently coded.
%  The parfor loop makes a big difference so it's worth keeping in this
%  form.
% This function works if both chains use smooth approximations to the ODE
% or if one or both of them uses the actual ODE solution.
% Use parallel_temper_tp if using parallel tempering in which case this
% is the divisor of the log top of the acceptance ratio.


% old_pars = {chainpars{chains_swap(1)}(:,iter),chainpars{chains_swap(2)}(:,iter)}
% old_X = {X{chains_swap(1)},X{chains_swap(2)}}
% log_r_bot = [log_r_bot{chains_swap(1)}(iter),log_r_bot{chains_swap(2)}(iter)]
% swap_counter=0

try
    
    for tt=1:2
        [log_r_top(tt) X_temp{tt}] = LogPostProb(Nknots,tspan,tgrid,prop_pars{tt},prop_X{tt},time,data,...
            prior_pars,odefn,tp(tt));
        %,old_xhvar{3-tt},old_D{3-tt},old_dvar{3-tt},old_sigsq{3-tt},old_xklogprob{3-tt});
    end
    
    if isinf(log_r_bot(2)) && ~isinf(log_r_bot(1))
        % if the chain with more flexibility (cooler temp) is okay but the higher
        % temp chain is stuck with isinf(log_r_bot)==1, then we should assume that
        % the lower temp chain is obligated to pass on its parameter value to the
        % higher temp chain.
        newpars=prop_pars;
        new_X = prop_X;
        swap_counter=swap_counter+1;
        log_r_bot=log_r_top;
    else
        if(exp(sum(log_r_top)-sum(log_r_bot))>rand )
            newpars=prop_pars;
            new_X = prop_X;
            swap_counter=swap_counter+1;
            log_r_bot=log_r_top;
        else
            newpars=old_pars;
            new_X = old_X;
        end
    end
    
catch
    newpars=old_pars;
    new_X = old_X;
    warning('problem in swap')
end

end

