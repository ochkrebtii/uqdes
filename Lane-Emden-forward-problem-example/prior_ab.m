function p_ab = prior_ab(ab,mu_ab,var_ab)
%   Setting up the prior for a,b as independent normals

p_ab = normpdf(ab,mu_ab,var_ab);
