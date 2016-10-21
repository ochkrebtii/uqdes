function [ Cov ] = RBF_II( hyperparameters, IntT1, IntT2 )


sigmaf = hyperparameters(1);
gamma  = hyperparameters(2); % Characteristic length parameter

lengthX = size(IntT1, 2);
lengthY = size(IntT2, 2);

%XY  = T(ones(lengthY,1),:)' - DerivT(ones(lengthX,1),:);
XY  = IntT2(ones(lengthX,1),:) - IntT1(ones(lengthY,1),:)';
XY2 = XY.^2;


%K = sigmaf*exp( (-1/(2*gamma^2)) * XY2 );

%Cov = -(1/gamma^2)*(XY).*K;

Cov = - ((sqrt(2)*gamma*sigmaf)/2).*( (XY./(sqrt(2)*gamma)) .* erf(XY./(sqrt(2)*gamma)) + ( -1 + exp(-(XY./(sqrt(2)*gamma)).^2) )./sqrt(pi) )...
      + ((sqrt(2)*gamma*sigmaf)/2)*( (IntT2(ones(lengthX,1),:)./(sqrt(2)*gamma)) .* erf(IntT2(ones(lengthX,1),:)./(sqrt(2)*gamma)) + ( -1 + exp(-(IntT2(ones(lengthX,1),:)./(sqrt(2)*gamma)).^2) )./sqrt(pi) )...
      + ((sqrt(2)*gamma*sigmaf)/2)*( (IntT1(ones(lengthY,1),:)'./(sqrt(2)*gamma)) .* erf(IntT1(ones(lengthY,1),:)'./(sqrt(2)*gamma)) + ( -1 + exp(-(IntT1(ones(lengthY,1),:)'./(sqrt(2)*gamma)).^2) )./sqrt(pi) );


%{
% Numerical
Cov_Num = zeros(lengthX, lengthY);
for i = 1:lengthX
    for j = 1:lengthY
        t1 = IntT1(i);
        t2 = IntT2(j);
        
        Cov_Num(i,j) = - ((sqrt(2)*gamma*sigmaf)/2)*( ((t2-t1)./(sqrt(2)*gamma)) * erf((t2-t1)./(sqrt(2)*gamma)) + ( -1 + exp(-((t2-t1)./(sqrt(2)*gamma)).^2) )./sqrt(pi) )...
                       + ((sqrt(2)*gamma*sigmaf)/2)*( (t2./(sqrt(2)*gamma)) * erf(t2./(sqrt(2)*gamma)) + ( -1 + exp(-(t2./(sqrt(2)*gamma)).^2) )./sqrt(pi) )...
                       + ((sqrt(2)*gamma*sigmaf)/2)*( (t1./(sqrt(2)*gamma)) * erf(t1./(sqrt(2)*gamma)) + ( -1 + exp(-(t1./(sqrt(2)*gamma)).^2) )./sqrt(pi) );
    end
end
%}


end