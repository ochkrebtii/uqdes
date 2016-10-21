function [logprobs X_prop um uv] = lik_ode(Nknots,tspan,tgrid,odefn,pars_prop,X_prop,time,data,um,uv)


DEpars          = pars_prop(1:2);
xfirstu         = pars_prop(3);
xfirstv         = pars_prop(4);
xlastu          = pars_prop(5);
alpha           = pars_prop(6);
lambda          = pars_prop(7);

x0              = [xfirstu,xfirstv];
nsolves         = 1;


if(nargin==8)
    [tgrid,um,uv] = podes_withvar(tspan,nsolves,Nknots,lambda,alpha,odefn,x0,DEpars,tgrid);
    cov = triu(uv)+triu(uv,1)';
    X_prop(:,1) = mvnrnd(um(:,1),cov);
    X_prop(:,2) = mvnrnd(um(:,2),cov);
else
    cov = triu(uv)+triu(uv,1)';
end


if ~isnan(sum(sum(X_prop)))
    dev              = data(1) - X_prop(end,1);
    sig              = sqrt(cov(end,end));
    try
        logprobs = - 0.5*(dev./sig).^2 - log(sqrt(2*pi) .* sig);
    catch
        logprobs = -Inf;
    end
else
    logprobs = -Inf;
end

