%% Demonstration of UQDES algorithm for the Lorenz system

clear all; close all; clc;

u0 = [-12,-5, 38]; sspan = [0 20]; theta = [10,8/3,28];
odefn = @lorenzode; nsolves = 1000; kernel = 'sqexp';  % sqexp or uniform
N = 2001; ds = range(sspan)/(N-1); lambda = 0.5*ds; alpha = 1;

disp('Please wait... ')
disp('(to speed up try, for e.g. uniform kernel with lambda=0.5ds)')

[uensemble,t,lnft] = uqdes(sspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta);
[ueuler,teuler] = euler(sspan,N,odefn,u0,theta);
[temp,trk4] = rk4(odefn,sspan(1),sspan(2),u0,N,theta);

disp(['UQDES algorithm run time: ', num2str(exp(lnft)/60),' minutes'])


%% Make Figure 2 (top and bottom panels)

stops = floor([0.35,0.485,0.56,1]*N);

figure(1)
stnames = {'u^(^1^)','u^(^2^)','u^(^3^)'};
for state = 1:length(u0);
    sp = subplot(length(u0),1,state);
    hold on
    for pp = 1:nsolves
        xflip = [t fliplr(t)];
        yflip = [uensemble(:,state,pp)' fliplr(uensemble(:,state,pp)')];
        p = patch(xflip,yflip,'r','EdgeAlpha',0.01,'FaceColor','none','linewidth',1.5);
    end
    axis tight
    if state == 3; xlabel('t'); end;
    ylabel(stnames{state})
    if state == 1
        ylims = get(sp,'YLim');
        for jj = 1:length(stops)
            line([t(stops(jj)),t(stops(jj))],ylims,'color','red','linewidth',2)
        end
    end
end


figure(2)
for kp = 1:4
    subaxis(1,4,kp) % can use subplot if you do not have subaxis package
    lp = stops(kp);
    plot3(reshape(uensemble(lp,1,:),nsolves,1),reshape(uensemble(lp ,2,:),nsolves,1),reshape(uensemble(lp ,3,:),nsolves,1),'.','markersize',8);
    view([45,140,60]);
    grid on
xlim([floor(min(reshape(uensemble(:,1,:),nsolves*N,1))),...
        ceil(max(reshape(uensemble(:,1,:),nsolves*N,1)))])
 ylim([floor(min(reshape(uensemble(:,2,:),nsolves*N,1))),...
        ceil(max(reshape(uensemble(:,2,:),nsolves*N,1)))])
 zlim([floor(min(reshape(uensemble(:,3,:),nsolves*N,1))),...
        ceil(max(reshape(uensemble(:,3,:),nsolves*N,1)))])
if kp == 3; xlabel('u^(^1^)'); end
if kp == 1; ylabel('u^(^2^)'); end
if kp == 1; zlabel('u^(^3^)'); end
title(['t = ',num2str(num2str(round(t(lp),2)))]);
end


figure(3)
for kp = 1:4
    subplot(1,4,kp)
    lp = stops(kp);
    [f,xi] = ksdensity(reshape(uensemble(lp,1,:),nsolves,1));
    %[yj,xj] = hist(reshape(xmean(lp,1,:),nsolves,1),50)
    plot(xi,f,'LineWidth',2)
    xlim([floor(min(reshape(uensemble(:,1,:),nsolves*N,1))),...
           ceil(max(reshape(uensemble(:,1,:),nsolves*N,1)))])
    xlabel(['u^(^1^)(',num2str(round(t(lp),2)),')']);
   if kp == 1; ylabel('relative frequency'); end
end

