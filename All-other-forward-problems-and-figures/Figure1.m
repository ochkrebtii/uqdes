%% Demonstration of uqdes algorithm

clear all; close all; clc;

u0 = [-1 0]; sspan = [0 10]; theta = 2;
odefn = @simpleode; odesoln = @simpleode_solution;
nsolves = 100; kernel = 'sqexp';  % sqexp or uniform
Nvec = [25,50,100];

figure
state = 2;
for nind = 1:length(Nvec)
    subaxis(1,length(Nvec),nind) % you may also use subplot here
    N = Nvec(nind);
    ds = range(sspan)/(N-1);
    lambda = 1*ds; alpha = N/100;
    [uensemble,t] = uqdes(sspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta);
    [ueuler,teuler] = euler(sspan,N,odefn,u0,theta);
    truth = odesoln(t,theta);
    truth = truth(1:2,:)';
    eu = plot(teuler,ueuler(:,state),'g--');
    hold on
    for pp = 1:nsolves
        tempt = sspan(1):0.01:sspan(2);
        tempensemble = spline(t,uensemble(:,state,pp)',tempt);
        xflip       = [tempt fliplr(tempt)];
        yflip       = [tempensemble fliplr(tempensemble)];
        p           = patch(xflip,yflip,'r','EdgeAlpha',0.025,'FaceColor','none','linewidth',1.5);
    end
    eu = plot(teuler,ueuler(:,state),'g--');
    tr = plot(t,truth(:,state),'r-');
    xlabel('t')
    if nind == 1; ylabel('u'); end
    axis([sspan(1),sspan(2),-4,4])
    box off
    %legend([eu,p,tr],'euler','uqdes','true')
end