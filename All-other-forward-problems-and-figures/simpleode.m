function r = simpleode(t,u,theta)

%% simpleode.m
%
% Second order ODE initial value problem with known closed form solution
% to be used for simulations
%
% Inputs
% t is time
% u is a 2x1 vector of states
% theta is a single parameter (e.g. try theta = 2)
% 
% Outputs
% r is the time derivative of u at time t
%

r    =  u;

if ndims(u) == 2

    r(1) = u(2);
    r(2) = -u(1)+sin(theta.*t);
    
elseif ndims(u) == 3
    
    r(:,1,:) = u(:,2,:)*ones(size(t));
    r(:,2,:) = -u(:,1,:)*ones(size(t))+sin(theta.*t);
    
end
