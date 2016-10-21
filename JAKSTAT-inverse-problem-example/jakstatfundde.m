function [r EpoRa] = jakstatfundde(t,y,Z,pars)%,randNum)

%% JakStat DDE in vector form

% if nargin == 4
%     randNum         = randn;
% end

% alphaepo            = pars.model(8);
% lambdaepo           = pars.model(9);

% if isinf(alphaepo)
    EpoRa           = interp1(pars.EpoRAtime,pars.EpoRAobs,t,'linear');
% else
%     [mean cov]      = interp(pars.EpoRAtime,pars.EpoRAobs,t,lambdaepo,alphaepo);
%     EpoRa           = mean + randNum*sqrt(cov);
% end

p                   = pars.model;

r                   =  y;

r(:,1)              = -p(1).*y(:,1).*EpoRa + 2*p(4).*Z(:,4);
r(:,2)              =  p(1).*y(:,1).*EpoRa - p(2).*y(:,2).^2;
r(:,3)              = -p(3).*y(:,3)        + 0.5.*p(2).*y(:,2).^2;
r(:,4)              =  p(3).*y(:,3)        - p(4).*Z(:,4);


end