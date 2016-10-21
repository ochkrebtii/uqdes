function r = lorenzode(t,u,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lorenzode
% The Lorenz system with three states 
% p is a 3x1 vector of parameters
% chaotic dynamics show up under, e.g. 
% p = [a,b,r] = [10,8/3,28];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r    = u;
r(1) = -p(1).*u(1) + p(1).*u(2);
r(2) = p(3).*u(1) - u(2) - u(1).*u(3);
r(3) = -p(2).*u(3) + u(1).*u(2);


end