function pj = prior_jeffrys(x)

% log uniform (Jeffrey's prior)
% Somethimes this is bad in which case an inverse Gamma is recommended and
% still conjugate.

pj = 1./x;
