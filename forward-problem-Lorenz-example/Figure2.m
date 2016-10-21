

% Currently set up with the same settings as the paper EXCEPT more draws


%% set up ODE model

% set random seed (or comment out to generate one automatically)
%RandStream.setGlobalStream(RandStream('mt19937ar','seed',83))

% Define ODE function and model parameters

%extreme3 = xmean(logical(sum(reshape(xmean(:,3,:)==max(max(xmean(:,3,:))),5000,100),2)),:,logical(sum(reshape(xmean(:,3,:)==max(max(xmean(:,3,:))),5000,100))))

inputs              = [10,8/3,28]';
xstart              = [-12,-5, 38];


odefn               = @lorenzode;

N                   = 3001;
nsolves             = 1000;
a                   = 0;
b                   = 20;
tspan               = [a b];

nstates             = 3;
nevalpts            = 1000;
evalgrid            = linspace(a,b,nevalpts);

lambda              = 2*(b-a)/(N-1);
alpha               = N;

disp('Go get yourself a cup of coffee, 1000 sampled solutions on a grid of 5001 points takes ~12 minutes ... ')

tic

[t,xmean,xcov]      = PODES_simple_fast(tspan,nsolves,N,lambda,alpha,odefn,xstart,inputs,evalgrid);

toc
%12.5 minutes




stops = [500,750,800,1000];

figure('Color',[1 1 1]);
stnames = {'u^(^1^)(t)','u^(^2^)(t)','u^(^3^)(t)'};
for state = 1:3;
    subaxis(3,1,state)
    hold on
    for kk = 1:nsolves
        xflip       = [t fliplr(t)];
        yflip       = [xmean(:,state,kk)' fliplr(xmean(:,state,kk)')];
        p           = patch(xflip,yflip,'r','EdgeAlpha',0.075,'FaceColor','none','Linesmoothing','on');
    end
    ylabel(stnames{state},'fontsize',20)
    if state == 1
        hold on
        for kp = 1:4
            plot(t(stops(kp))*[1,1],[-20,20],'-r','LineWidth',3)
        end
        hold off
    end
    box on
    if state == 1
        axis([a,b,-20,20])
    elseif state == 2
        axis([a,b,-30,30])
    else
        axis([a,b,0,50])
    end
    set(gca,'Linewidth',1.5);
    hold off
end

xlabel('t','fontsize',20)
%xlabel(strcat(num2str(lambda),'alpha',num2str(alpha)),'fontsize',20)

%
figure
for kp = 1:4
    subplot(1,4,kp)
    lp = stops(kp);
    [f,xi] = ksdensity(reshape(xmean(lp,1,:),nsolves,1))
    %[yj,xj] = hist(reshape(xmean(lp,1,:),nsolves,1),50)
    plot(xi,f,'LineWidth',3)
    xlim([floor(min(reshape(xmean(:,1,:),nsolves*nevalpts,1))),
           ceil(max(reshape(xmean(:,1,:),nsolves*nevalpts,1)))])
    title(strcat(' \fontsize{16} Time = ',num2str(t(lp))))
end


figure
for(kp = 1:4)
    subplot(1,4,kp)
    lp = stops(kp)
plot3(reshape(xmean(lp ,1,:),nsolves,1),reshape(xmean(lp ,2,:),nsolves,1),reshape(xmean(lp ,3,:),nsolves,1),'.');grid on
 xlim([floor(min(reshape(xmean(:,1,:),nsolves*nevalpts,1))),
        ceil(max(reshape(xmean(:,1,:),nsolves*nevalpts,1)))])
 ylim([floor(min(reshape(xmean(:,2,:),nsolves*nevalpts,1))),
        ceil(max(reshape(xmean(:,2,:),nsolves*nevalpts,1)))])
 zlim([floor(min(reshape(xmean(:,3,:),nsolves*nevalpts,1))),
        ceil(max(reshape(xmean(:,3,:),nsolves*nevalpts,1)))])

xlabel('\fontsize{14}u^(^1^)(t)');
ylabel('\fontsize{14}u^(^2^)(t)');
zlabel('\fontsize{14}u^(^3^)(t)')
title(strcat(' \fontsize{16} Time = ',num2str(t(lp))) )
end




