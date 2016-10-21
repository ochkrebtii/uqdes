function pri=prior_chisq(x,df)
% prior for the initial values.  


pri = chi2pdf(x,df);
