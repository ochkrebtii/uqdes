function p_ab = prior_normal(x,mu,var)
%   Setting up the prior for a,b as independent normals

p_ab = normpdf(x,mu,sqrt(var));
