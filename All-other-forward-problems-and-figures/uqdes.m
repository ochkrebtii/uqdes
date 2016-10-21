function [uensemble,t,logfntime,m_deriv_svec] ...
    = uqdes(sspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta,tgrid)

%% uqdes.m 
%
% DESCRIPTION
% Bayesian updating algorithm by Oksana A. Chkrebtii last updated
% May 20, 2016. The probabilistic solver is described int he paper Bayesian
% Uncertainty Quantification for Differential Equations, by O.A. Chkrebtii,
% D.A. Campbell, B. Calderhead, M.A. Girolami. This implementation allows
% the option to evaluate the forward model outside of the discretization
% grid.
%
% INPUTS
% sspan = 1x2 vector containing the upper and lower bound of the domain of
% integration
% N = an integer for the discretizatio mesh size
% odefn = a function handle for the system of ODEs
% u0 = 1xM vector of initial conditions
% theta = parameters used in the system of ODEs, dimension of theta depends
% on the inputs of the function handle odefn
% tgrid = 1xq optional parameter passed if we wish to evaluate the forward
% model at time points lying outside the discretization grid
% 
% OUTPUTS
% This function returns 
% t: the discretization grid at which the function was evaluated (for the
% moment these are equally spaced)
% logftime: gives the log computation time in seconds
% 
% EXAMPLE
% [umat,tvec,lft] = uqdes([0,10],5,50,'sqexp',1,50,@simpleode,[-1,0],2);


tic % start timing  

% preallocate variables
M = length(u0); B = nsolves; 
s = linspace(sspan(1),sspan(2),N); ds = s(2)-s(1);
if nargin == 10 
    t = sort(unique([s,tgrid]));
    [~,tinds] = intersect(t,s);
    [~,sinds] = intersect(t,tgrid);
else t = s; tinds = (1:N)'; sinds = []; end; 
if size(u0,1)>size(u0,2); u0 = u0'; end;

switch lower(kernel) % sort out which kernel to use
    case {'sqexp'};  kern = 'se'; trim = N;
    case {'uniform'}; kern = 'un'; trim = 2*ceil(lambda/ds);
    otherwise; disp('Unknown kernel -- try agaian'); return;
end

% define appropriate kernel convolutions
QQ = str2func(strcat('QQ1d_',kern));
RR = str2func(strcat('RR1d_',kern));
QR = str2func(strcat('QR1d_',kern));

uensemble = repmat(u0,[length(t),1,B]);
f = odefn(s(1),uensemble(1,:,:),theta);
% there are two possibilities for the prior mean: zero derivative, or
% constant exact derivative f(s(1),u0): either one works
m_deriv_svec = repmat(0*u0,[length(t),1,B]);
%m_deriv_svec = repmat(f,[length(t),1,1]);
m_state_svec =  uensemble + bsxfun(@times,m_deriv_svec,repmat(t',1,M));
C_deriv_ssmat = RR(t,t,lambda,sspan(1),sspan(2))/alpha;
C_state_ssmat = QQ(t,t,lambda,sspan(1),sspan(2))/alpha;
C_cross1_ssmat = QR(t,t,lambda,sspan(1),sspan(2))/alpha;
kinv = 1/(C_deriv_ssmat(1,1));
f_diff = kinv*(f-m_deriv_svec(1,:,:));      
randnNums  = randn(length(t),M,B);  % generate random numbers outside of loop
counter = 0;

% run one-step uqdes algorithm
for n = tinds(1:end-1)'
    ind = tinds(max(counter+1,max(1,counter-trim)):min(end,min(end,counter+trim)));
	endind = vertcat(tinds(counter+1:N),sinds);
    counter = counter + 1;
    nextind = tinds(counter+1);
    m_state_svec(endind,:,:) = m_state_svec(endind,:,:) + bsxfun(@times,C_cross1_ssmat(endind,n),f_diff(1,:,:));
    m_deriv_svec(ind,:,:) = m_deriv_svec(ind,:,:) + bsxfun(@times,C_deriv_ssmat(ind,n),f_diff(1,:,:));
    C_state_ssmat(endind,endind) = C_state_ssmat(endind,endind) - kinv*C_cross1_ssmat(endind,n)*(C_cross1_ssmat(endind,n))';
    C_cross1_ssmat(endind,ind) = C_cross1_ssmat(endind,ind) - kinv*C_cross1_ssmat(endind,n)*C_deriv_ssmat(n,ind);
    C_deriv_ssmat(ind,ind) = C_deriv_ssmat(ind,ind) - kinv*C_deriv_ssmat(ind,n)*C_deriv_ssmat(n,ind);
    uensemble(nextind,:,:) = m_state_svec(nextind,:,:) + sqrt(C_state_ssmat(nextind,nextind))*randnNums(n,:,:);   
    kinv = 1/(C_deriv_ssmat(nextind,nextind)+C_deriv_ssmat(tinds(counter),tinds(counter)));
    f_diff = kinv*(odefn(t(nextind),uensemble(nextind,:,:),theta) - m_deriv_svec(nextind,:,:));
    if sum(f_diff >= 1e10 | isnan(f_diff) | isinf(f_diff))>0
        disp('Algorithm failed to converge: try increasing mesh size or changing assumptions')
        uensemble = [];
        logfntime = [];
        m_deriv_svec = [];
        return
    end
end
    
if nargin == 10
    randnNums  = randn(length(sinds),M,B);  % generate random numbers outside of loop
    uensemble(sinds,:,:) = m_state_svec(sinds,:,:) + bsxfun(@times,randnNums,sqrt(diag(C_state_ssmat(sinds,sinds))));   
end

logfntime = log(toc); % end timer


