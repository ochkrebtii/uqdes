%% Demonstration of uqdes algorithm


clear all; close all; clc;
u0 = [-1 0]; sspan = [0 10]; theta = 2; 
odefn = @simpleode; odesoln = @simpleode_solution;
nsolves = 10; kernel = 'uniform';  % choose between sqexp or uniform
Nvec = floor(linspace(50,1000,5)); 
lambdafactorvec = linspace(1,5,3); 
logftuqdes = zeros(length(lambdafactorvec),length(Nvec));
logimsevecuqdes = zeros(length(lambdafactorvec),length(Nvec));
logfteuler = zeros(1,length(Nvec));
logimseveceuler = zeros(1,length(Nvec));

for nn = 1:length(Nvec)
    
    try
    nevalpts = Nvec(nn);
    N = Nvec(nn);
    ds = range(sspan)/(N-1);
    lambdavec = lambdafactorvec*ds;
    alpha = N; % scale the prior precision with respect to the number of grid poitns used
    t = linspace(sspan(1),sspan(2),N);
    truth = odesoln(t,theta);
    truth = truth(1:2,:)';
    
    for ll = 1:length(lambdavec)
        lambda = lambdavec(ll);
        [ueuler,teuler,logfteuler(ll,nn)] = euler(sspan,N,odefn,u0,theta);
        [uuqdes,tuqdes,logftuqdes(ll,nn)] = uqdes(sspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta);
        logimsevecuqdes(ll,nn) = log(imse(tuqdes,uuqdes(:,:,:),truth));
        logimseveceuler(ll,nn) = log(imse(teuler,ueuler,truth));
    end
    
    catch
    end

end

%% Plot Figure

figure
s2 = subplot(1,2,1);
for ll=1:length(lambdavec)
    h(ll) = plot(Nvec,logimsevecuqdes(ll,:),'linestyle','-');
    hold on
    legInfo{ll} = ['\lambda/h = ',num2str(lambdavec(ll)/ds)];
end
h(ll+1) = plot(Nvec,logimseveceuler(1,:),'r--');
legInfo{ll+1} = 'Euler';
xlabel('N'); ylabel('log IMSE');
axis tight; box on;
legend(h,legInfo,'Location','NorthEast','Orientation','Horizontal')
s1 = subplot(1,2,2);
plot(Nvec,mean(logftuqdes,1),'linestyle','-')
hold on
plot(Nvec,mean(logfteuler,1),'linestyle','--')
xlabel('N'); ylabel('log time');
axis tight
legend('UQDES-1 (U) 10 draws','Euler-1', 'Location','SouthEast','Orientation','Vertical')

