function pri=prior_exp(x,lam)

% exponential prior with mean parameter mu = 1/lam

pri = exppdf(x,1./lam);
