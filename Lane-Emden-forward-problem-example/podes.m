function [u_ensemble,evalgrid,function_time] = podes(tspan,nsolves,N,kernel,lambda,alpha,odefn,u0,theta,evalgrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PODES implementation by Oksana A. Chkrebtii last updated on April 27,
% 2014. The probabilistic solver is described int he paper Bayesian
% Uncertainty Quantification for Differential Equations, by O.A. Chkrebtii,
% D.A. Campbell, M.A. Girolami, B. Calderhead.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic             % time this function

M               = length(u0);
B               = nsolves;
a               = tspan(1);
b               = tspan(2);
t               = linspace(a,b,N);

if size(u0,1) > size(u0,2)
    u0 = u0';
end

if nargin == 9
    evalgrid = t;
end

u               = repmat(vertcat(u0,zeros(N-1,M)),[1,1,B]);
Du              = repmat(vertcat(odefn(t(1),u0,theta),zeros(N-1,M)),[1,1,B]);
u_ensemble      = repmat(repmat(u0,length(evalgrid),1),[1,1,B]);


switch lower(kernel)
   case {'gaussian'}
      disp('Kernel: square exponential')
      disp('(note you have assumed an infinitely differentiable smooth solution)')
      kern        = 'se';
      Ntrim       = N;
   case {'uniform'}
      disp('Kernel: uniform')
      disp('(note that you have assumed a possibly discontinuous first derivative)')
      kern        = 'un';
      Ntrim       = ceil(2*lambda*N/range(tspan))+1;
   otherwise
      disp('Unknown kernel!!!');
      return
end

% define appropriate kernel convolutions
QQ = str2func(strcat('QQ1d_',kern));
RR = str2func(strcat('RR1d_',kern));
RQ = str2func(strcat('RQ1d_',kern));
QR = str2func(strcat('QR1d_',kern));


% Generate random numbers outside of loop
randnNums                    = randn(B,M,N);


for k = 1:N
    
%     if mod(k, 1000) == 0
%         disp(k)
%     end
    
    if k == 1
        Binv(1,1)              = alpha/RR(t(1),t(1),lambda,a,b);
    else
        % components of the inverse of block matrix (Sigma + RR)^{-1}
        if k < Ntrim + 1
            b                       = RR(t(1:k-1),t(k),lambda,a,b)./alpha;
            D                       = Binv(1:k-1,1:k-1);
            Db                      = D*b;
            btD                     = Db';
            c                       = ss_kp1 + RR(t(k),t(k),lambda,a,b)./alpha - b'*Db;
        else
            bt                      = RR(t(Ntrim),t(1:Ntrim-1),lambda,a,b)./alpha;
            D                       = Binv(k-Ntrim+1:k-1,1:k-1);
            btD                     = bt*D;
            Db                      = btD';
            c                       = ss_kp1 + RR(t(k),t(k),lambda,a,b)./alpha - btD(1,k-Ntrim+1:k-1)*bt';
        end
        
        % build up the inverse of block matrix (Sigma + RR)^{-1}
        Binv                      = Binv + Db*(btD./c);
        Binv(1:k-1,k)             = -Db./c;
        Binv(k,1:k-1)             = -btD./c;
        Binv(k,k)                 = 1/c;
    end
    
    if k<N
        
        qrb                          = ( QR(t(k+1),t(1:k),lambda,a,b)*Binv )./alpha;
        ss_kp1                       = RR(t(k+1),t(k+1),lambda,a,b)./alpha - (RR(t(k+1),t(1:k),lambda,a,b)*Binv*RR(t(1:k),t(k+1),lambda,a,b) )/(alpha^2);
        sapv                         = QQ(t(k+1),t(k+1),lambda,a,b)/alpha - ( qrb*RQ(t(1:k),t(k+1),lambda,a,b) )./alpha;
           
        for bb = 1:B
            sapm                     = u(1,:,bb) + qrb*Du(1:k,:,bb);
            u(k+1,:,bb)              = sapm + randnNums(bb,:,k)*sqrt(sapv);
            Du(k+1,:,bb)             = odefn(t(k+1),u(k+1,:,bb),theta);
        end
    else
        
        qrbeval                      = QR(evalgrid,t(1:k),lambda,a,b)*(Binv./alpha);
        
        for bb = 1:B
            u_ensemble(:,:,bb)       = repmat(u(1,:,bb),length(evalgrid),1) + qrbeval*Du(1:k,:,bb);
        end
        
    end
   
end

function_time = toc;    % time this function

end