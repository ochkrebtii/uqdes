function r = lorenzode(t,u,p)

%% lorenzode.m
%
% The Lorenz system with three states
%
% Inputs
% t is time
% u is a 3x1 vector of states
% p is a 3x1 vector of parameters (e.g. p = [10,8/3,28]' for chaotic)
% 
% Outputs
% r is the 3x1 gradient of u
%

r = u; 

if ndims(u) == 2
    
    r(1) = -p(1).*u(1) + p(1).*u(2);
    r(2) = p(3).*u(1) - u(2) - u(1).*u(3);
    r(3) = -p(2).*u(3) + u(1).*u(2);
    
    
elseif ndims(u) == 3
    
    r(:,1,:) = (-p(1).*u(:,1,:) + p(1).*u(:,2,:))*ones(size(t));
    r(:,2,:) = (p(3).*u(:,1,:) - u(:,2,:) - u(:,1,:).*u(:,3,:))*ones(size(t));
    r(:,3,:) = (-p(2).*u(:,3,:) + u(:,1,:).*u(:,2,:))*ones(size(t));

end
