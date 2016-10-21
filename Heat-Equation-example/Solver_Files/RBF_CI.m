function [ Cov ] = RBF_CI( hyperparameters, T, IntT )


sigmaf     = hyperparameters(1);
charlength = hyperparameters(2);

lengthX = size(T, 2);
lengthY = size(IntT, 2);

%XY  = T(ones(lengthY,1),:)' - DerivT(ones(lengthX,1),:);
XY  = IntT(ones(lengthX,1),:) - T(ones(lengthY,1),:)';
XY2 = XY.^2;


%K = sigmaf*exp( (-1/(2*charlength^2)) * XY2 );

%Cov = -(1/charlength^2)*(XY).*K;

Cov = (sigmaf/2)*(erf(XY./(sqrt(2)*charlength)) + erf(T(ones(lengthY,1),:)'./(sqrt(2)*charlength)) );


end