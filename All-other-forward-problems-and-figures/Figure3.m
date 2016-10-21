%% Demonstration of uqdes algorithm

clear all; close all;

u0 = [-1 0]; sspan = [0 10]; theta = 2;
odefn = @simpleode; odesoln = @simpleode_solution;
nsolves = 5; N = 50; ds = range(sspan)/(N-1);
lambda = 1*ds; alpha = N/10; kernel = 'sqexp';  % sqexp or uniform

[uensemble,t,function_time] ...
    = plot_uqdes(sspan,N,kernel,lambda,alpha,odefn,odesoln,u0,theta);
 

% %% verify that the ensemble looks reasonable
% figure
% truth = odesoln(t,theta);
% truth = truth(1:2,:)';
% subplot(1,2,1)
% plot(t,euler(:,1),'g-')
% hold on
% for pp = 1:nsolves
%     plot(t,uensemble(:,1,pp),'sb')
% end
% plot(t,truth(:,1),'r--')
% legend('euler','uqdes','true')
% subplot(1,2,2)
% plot(t,euler(:,2),'g-')
% hold on
% for pp = 1:nsolves
%     plot(t,uensemble(:,2,pp),'sb')
% end
% plot(t,truth(:,2),'r--')
% legend('euler','uqdes','true')
