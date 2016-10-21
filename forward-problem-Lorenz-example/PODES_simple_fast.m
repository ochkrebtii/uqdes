function [evalgrid,um,uv] = PODES_simple_fast(tspan,nsolves,Nknots,lambda,alpha,odefn,x0,pars,evalgrid)


N               = Nknots;
M               = length(x0);
B               = nsolves;
a               = tspan(1);
b               = tspan(2);
t               = linspace(a,b,N);

if nargin == 8
    evalgrid = t;
end

u               = repmat(vertcat(x0,zeros(N-1,M)),[1,1,B]);
Du              = repmat(vertcat(odefn(t(1),x0,pars),zeros(N-1,M)),[1,1,B]);
um              = repmat(repmat(x0,length(evalgrid),1),[1,1,B]);
meaneval        = um;

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
        ss_kp1                      = RR1d(t(k+1),t(k+1),lambda,a,b)./alpha - (RR1d(t(k+1),t(1:k),lambda,a,b)*Binv*RR1d(t(1:k),t(k+1),lambda,a,b) )/(alpha^2);
        sapv                        = QQ1d(t(k+1),t(k+1),lambda,a,b)/alpha - ( qrb*RQ1d(t(1:k),t(k+1),lambda,a,b) )./alpha;
   
        randnNums                    = randn(B,M); % Generate random numbers outside of loop
        
        for bb = 1:B
            sapm                     = u(1,:,bb) + qrb*Du(1:k,:,bb);
            u(k+1,:,bb)              = sapm + randnNums(bb,M)*sqrt(sapv);
            Du(k+1,:,bb)             = odefn(t(k+1),u(k+1,:,bb),pars);
        end
    else
        
        qrbeval                     = QR1d(evalgrid,t(1:k),lambda,a,b)*(Binv./alpha);
        uv                          = QQ1d(evalgrid,evalgrid,lambda,a,b)./alpha - qrbeval*(RQ1d(t(1:k),evalgrid,lambda,a,b)./alpha);
        
        for bb = 1:B
            um(:,:,bb)              = repmat(u(1,:,bb),length(evalgrid),1) + qrbeval*Du(1:k,:,bb);
        end
        
    end
   
end


end