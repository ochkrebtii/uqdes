function r = decayivp(t,u,theta)

%% decayivp.m
%
% Decay ODE initial value problem with one state and known analytical
% solution, to be used for numerical experiments
% 
% Inputs
% t is time
% u is the state
% theta is a single parameter
% 
% Outputs
% r is the time derivative of u at time t
%

r    =  u;

if ndims(u) == 2
    
    r(1) = u(1).*(1-theta.*t);
    
elseif ndims(u) == 3
    
    r(:,1,:) = u(:,1,:).*(1-theta.*t);
    
end