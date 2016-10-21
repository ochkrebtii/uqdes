function p_al = prior_lognorm(x,shift,mu,var)

% Shifted log normal prior

p_al = lognpdf(x-shift,mu,sqrt(var));