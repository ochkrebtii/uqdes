function [sol] = simpleode_solution(t,theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact solution to our toy second order ODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = reshape(t,1,length(t));

solstate = (1/(theta^2-1))*(-(theta^2).*cos(t) + theta*sin(t) - sin(theta*t) + cos(t));
solderiv = (1/(theta^2-1))*((theta^2).*sin(t) + theta*cos(t) - theta*cos(theta*t) - sin(t));
solsecderiv = (1/(theta^2-1))*((theta^2).*cos(t) - theta*sin(t) + theta^2*sin(theta*t) - cos(t));

sol = vertcat(solstate,solderiv,solsecderiv);

end