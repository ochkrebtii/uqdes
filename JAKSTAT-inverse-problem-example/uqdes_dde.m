function [um,uv,evalgrid]...
    = uqdes_dde(tspan,nsolves,N,laglength,...
    lambda,alpha,ddefn,ddehist,pars,histpars,...
    evalgrid)

a               = tspan(1);
b               = tspan(2);
u0              = ddehist(a,histpars);
M               = length(u0);
B               = nsolves;
tau             = laglength;

t               = linspace(a,b,N);

if nargin < 11
    evalgrid    = t;
end

t = [t,0,0];

[histmean histsv histsd]        = ddehist(t(1)-tau,histpars);
u_lagged                        = repmat(histmean,B,1) + randn(B,M).*repmat(sqrt(histsv),B,1);   %obtain first lagged state

u                               = repmat(vertcat(u0,zeros(N-1,M)),[1,1,B]);
um                              = repmat(repmat(u0,length(evalgrid),1),[1,1,B]);
Du                              = zeros(N,M,B);

%randEpo                        = randn(B,1);
randnNums                       = randn(B,M,N); % Generate random numbers outside of loop
%randnNumsEpo                   = randn(B,M);

for bb = 1:B

    Du(1,:,bb) = ddefn(t(1),u0,u_lagged(bb,:),pars);%,randEpo(bb));

end

kernel          = 'gaussian';

if strcmp(kernel,'uniform')

    Ntrim       = ceil(2*lambda*N/range(tspan))+1;

else
    
    Ntrim       = N;

end


for k = 1:N
    
    if mod(k, 1000) == 0
        disp(k)
    end
    
    if k == 1
        Binv(1,1)              = alpha/RR1d(t(1),t(1),lambda,a,b);
    else
        % components of the inverse of block matrix (Sigma + RR1d)^{-1}
        if k < Ntrim + 1
            
            b                       = RR1d(t(1:k-1),t(k),lambda,a,b)./alpha;
            D                       = Binv(1:k-1,1:k-1);
            Db                      = D*b;
            btD                     = Db';
            c                       = ss_kp1 + RR1d(t(k),t(k),lambda,a,b)./alpha - b'*Db;
        else
            
            bt                      = RR1d(t(Ntrim),t(1:Ntrim-1),lambda,a,b)./alpha;
            D                       = Binv(k-Ntrim+1:k-1,1:k-1);
            btD                     = bt*D;
            Db                      = btD';
            c                       = ss_kp1 + RR1d(t(k),t(k),lambda,a,b)./alpha - btD(1,k-Ntrim+1:k-1)*bt';
            
        end
        
        % build up the inverse of block matrix (Sigma + RR1d)^{-1}
        Binv                      = Binv + Db*(btD./c);
        Binv(1:k-1,k)             = -Db./c;
        Binv(k,1:k-1)             = -btD./c;
        Binv(k,k)                 = 1/c;
        
    end
    
    if k<N
        
        qrb                         = ( QR1d(t(k+1),t(1:k),lambda,a,b)*Binv )./alpha;
        rrb_lag                     = ( RR1d(t(k+2)-tau,t(1:k),lambda,a,b)*Binv )./alpha;
        qrb_lag                     = ( QR1d(t(k+2)-tau,t(1:k),lambda,a,b)*Binv )./alpha;
        ss_kp1                      = RR1d(t(k+1),t(k+1),lambda,a,b)./alpha - (RR1d(t(k+1),t(1:k),lambda,a,b)*Binv*RR1d(t(1:k),t(k+1),lambda,a,b) )/(alpha^2)...
                                        + RR1d(t(k+2)-tau,t(k+2)-tau,lambda,a,b)./alpha - rrb_lag*RR1d(t(1:k),t(k+2)-tau,lambda,a,b)./alpha...
                                        + 2*RR1d(t(k+2)-tau,t(k+1),lambda,a,b)./alpha - rrb_lag*RR1d(t(1:k),t(k+1),lambda,a,b)/alpha;   % step-ahead pointwise variance for the derivative
        sapv                        = QQ1d(t(k+1),t(k+1),lambda,a,b)./alpha - ( qrb*RQ1d(t(1:k),t(k+1),lambda,a,b) )./alpha;            % step-ahead pointwise variance for the state
          
        
        
        for bb = 1:B
            
            sapm                     = u(1,:,bb) + qrb*Du(1:k,:,bb);
            tester = sum(sum(sapm));
            if isnan(tester) || isinf(tester)
            
                um = NaN;
                uv = NaN;
                return
            
            end
            
            u(k+1,:,bb)                 = sapm + randnNums(bb,:,k)*sqrt(sapv);
            Du(k+1,:,bb)                = ddefn(t(k+1),u(k+1,:,bb),u_lagged(bb,:),pars);%,randnNumsEpo(bb,M));
            
            if t(k+2)-tau > a
                
                u_lagged(bb,:)           = u(1,:,bb) + qrb_lag*Du(1:k,:,bb);
            
            else
                
                hmean                    = ddehist(t(k+2)-tau,histpars);
                u_lagged(bb,:)           = hmean; 
            
            end
            
        end
        
    else
        
        qrbeval                     = QR1d(evalgrid,t(1:k),lambda,a,b)*(Binv./alpha);
        uv                          = QQ1d(evalgrid,evalgrid,lambda,a,b)./alpha - qrbeval*(RQ1d(t(1:k),evalgrid,lambda,a,b)./alpha);
        
        for bb = 1:B
            
            um(:,:,bb)              = repmat(u(1,:,bb),length(evalgrid),1) + qrbeval*Du(1:k,:,bb);
            
        end
        
    end
end


