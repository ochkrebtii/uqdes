% Demonstration of uqdes algorithm
% Modify simulation settings below and run

%% Simulation settings

%s = RandStream('mt19937ar','Seed',33);
%RandStream.setGlobalStream(s);

%clear all; close all; clc;

u0 = [-1,0]; sspan = [0 10]; 
odefn = @simpleode; odesoln = @simpleode_solution;
N = 50; ds = range(sspan)/(N-1); t = linspace(sspan(1),sspan(2),N);
theta = 2; truth = odesoln(t,theta);  truth = truth(1:2,:)';
nsolves = 100; kernel = 'sqexp';  % choose between sqexp or uniform

lambdafactorvec = 1; lambda = lambdafactorvec*ds; alpha = N/10; 
[uuqdes,tuqdes,lftuqdes,duqdes] ...
    = uqdes(sspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta);
[ueuler,teuler,lfteuler,deuler] = euler(sspan,N,odefn,u0,theta);

disp(['uqdes IMSE: ', num2str(log(imse(tuqdes,uuqdes,truth)))])
disp(['euler IMSE: ', num2str(log(imse(teuler,ueuler,truth)))])


figure

subplot(1,2,1)
plot(t,truth(:,1),'r--')
hold on
plot(teuler,ueuler(:,1),'g--')
for pp = 1:nsolves
    plot(tuqdes,uuqdes(:,1,pp),'.b')
end
legend('true','euler','uqdes')
subplot(1,2,2)
plot(t,truth(:,2),'r--')
hold on
plot(teuler,ueuler(:,2),'g--')
for pp = 1:nsolves
    plot(tuqdes,uuqdes(:,2,pp),'.b')
end
legend('true','euler','uqdes')
