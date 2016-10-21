function [modpars,tau,histvec,alpha,lambda] = extractparms(pars)

%% This function extracts model and auxiliary parameters from the pars array

modpars              = pars(1:6);        % model parameters
tau                  = pars(7);          % legnth of time lag
histvec              = [pars(8);0;0;0];  % history
alpha                = pars(9);
lambda               = pars(10);
%alphaepo             = pars(11);
%lambdaepo            = pars(12);

%EpoRA               = horzcat(pars.EpoRAtime,pars.EpoRAobs);

end