function [ueuler,s,logfntime,deuler] = euler(sspan,N,odefn,u0,theta)

%% euler.m
% This function implements the one-step Euler method for approximating the
% solution of the differential equation system odefn
%
% Inputs
% sspan = 1x2 vector containing the upper and lower bound of the domain of
% integration
% N = an integer for the discretizatio mesh size
% odefn = a function handle for the system of ODEs
% u0 = 1xM vector of initial conditions
% theta = parameters used in the system of ODEs, dimension of theta depends
% on the inputs of the function handle odefn
% 
% Outputs
% This function returns 
% ueuler: a NxM vector of approximate ODE solution evaluated at the N grid
% points
% t: the discretization grid at which the function was evaluated (for the
% moment these are equally spaced)
% logftime: gives the log computation time in seconds
% 
% e.g. [umat,tvec,lft] = euler([0,10],50,[-1,0],odefn,2);


%% Set up some variables
M = length(u0);
s = linspace(sspan(1),sspan(2),N);
ds = s(2)-s(1);
if size(u0,1)>size(u0,2); u0 = u0'; end

%% Run One-step Euler method
tic % start timer
ueuler = repmat(u0,N,1);
deuler = repmat(odefn(s(1),u0,theta),N,1);
for n = 1:N-1
    deuler(n+1,:) = odefn(s(n+1),ueuler(n,:)',theta);
    for m = 1:M
        ueuler(n+1,m) = ueuler(n,m) + deuler(n+1,m)*ds;
    end
    ueuler(n+1:end,:) = repmat(ueuler(n+1,:),N-n,1);
end
logfntime = log(toc); % end timer
end

