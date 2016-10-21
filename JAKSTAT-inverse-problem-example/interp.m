function [mu cov] = interp(tobs,yobs,tpred,lambda,alpha)


yobs                            = reshape(yobs,length(yobs),1);
d                               = length(tobs);
nugget                          = 1e-10;

% for matrix inversion
L                               = chol( eye(d,d)*nugget + RR1d(tobs,tobs,lambda)/alpha );
A                               = L\(L'\yobs);

V                               = L'\RR1d(tobs,tpred,lambda)/alpha;

mu                              = RR1d(tpred,tobs,lambda)*A/alpha;
cov                             = RR1d(tpred,tpred,lambda)/alpha - V'*V;

end