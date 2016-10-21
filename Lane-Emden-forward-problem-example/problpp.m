function [lpp X_prop um uv] = problpp(Nknots,tspan,tgrid,pars_prop,X_prop,time,data,prior_pars,odefn,tp,um,uv)

%pars_prop = pars_start;
%X_prop=x_start;

xfirstu         = pars_prop(3);

if (nargin==10)
    [logprobs X_prop um uv] = lik_ode(Nknots,tspan,tgrid,odefn,pars_prop,X_prop,time,data);
else
    [logprobs X_prop um uv] = lik_ode(Nknots,tspan,tgrid,odefn,pars_prop,X_prop,time,data,um,uv);
end


if(~isnan(sum(sum(X_prop))))
    lpp = sum(tp*logprobs) + sum(log(prior_normal(xfirstu,prior_pars(1),prior_pars(2))));
    
    if ~isreal(lpp)
        lpp = -Inf;
    end
else
    lpp = -Inf;
end

end