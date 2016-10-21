% Demonstration of uqdes algorithm
% Modify simulation settings below and run

%% Simulation settings

%s = RandStream('mt19937ar','Seed',33);
%RandStream.setGlobalStream(s);

close all; clear all; clc;

u0 = 1; sspan = [0 5]; 
odefn = @decayivp; odesoln = @decayivp_solution; odederiv = @decayivp_derivative;
N = 100; ds = range(sspan)/(N-1); t = linspace(sspan(1),sspan(2),N);
theta = 2; truth = odesoln(t,theta); truth_deriv = odederiv(t,theta);
nsolves = 100; kernel = 'sqexp';  % choose between sqexp or uniform

lambdafactorvec = 0.5; lambda = lambdafactorvec*ds; alpha = N;
[uuqdes,tuqdes,lftuqdes,duqdes] ...
    = uqdes(sspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta);
[ueuler,teuler,lfteuler,deuler] = euler(sspan,N,odefn,u0,theta);
disp(['uqdes IMSE: ', num2str(log(imse(tuqdes,uuqdes,truth)))])
disp(['uqdes derivative IMSE: ', num2str(log(imse(tuqdes,duqdes,truth_deriv)))])
disp(['euler IMSE: ', num2str(log(imse(teuler,ueuler,truth)))])
disp(['euler derivative IMSE: ', num2str(log(imse(teuler,deuler,truth_deriv)))])

figure

subplot(1,2,1)
leg(1) = plot(teuler,ueuler(:,1),'g-');
hold on
for bb = 1:nsolves
    leg(2) = plot(tuqdes,uuqdes(:,1,bb),'.c');
end
leg(3) = plot(t,truth(:,1),'r--');
legend('euler','uqdes','true')
title('State')
axis tight

subplot(1,2,2)
pe = plot(teuler,deuler(:,1),'g-');
hold on
for bb = 1:nsolves
    pp = plot(tuqdes,duqdes(:,1,bb),'sb');
end
pt = plot(t,truth_deriv,'r--');
leginfo={'euler','uqdes','true'};
legend(leg,leginfo)
title('Derivative')
axis tight
