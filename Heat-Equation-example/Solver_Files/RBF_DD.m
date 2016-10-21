function [ Cov ] = RBF_DD( hyperparameters, DerivT1, DerivT2 )


sigmaf     = hyperparameters(1);
charlength = hyperparameters(2);

lengthX = size(DerivT1, 2);
lengthY = size(DerivT2, 2);

XY  = DerivT2(ones(lengthX,1),:) - DerivT1(ones(lengthY,1),:)';
XY2 = XY.^2;

%{
for a = 1:size(X, 2)
    for b = 1:size(Y, 2)
        K(a,b) = sigmaf*exp( (-1/(2*charlength^2)) * (X(a)-Y(b))^2 );
    end
end
%}
K = (sigmaf./(sqrt(2*pi)*charlength))*exp( (-1/(2*charlength^2)) * XY2 );

%{
for a = 1:size(X, 2)
    for b = 1:size(Y, 2)
        Kstar1(a,b) = -(sigmaf/charlength^2)*(X(a)-Y(b))*exp( (-1/(2*charlength^2)) * (X(a)-Y(b))^2 );
    end
end
%}
%Cov = (1/charlength^2)*(XY).*K;

%{
for a = 1:size(X, 2)
    for b = 1:size(Y, 2)
        starK(a,b) = (sigmaf/charlength^2)*(X(a)-Y(b))*exp( (-1/(2*charlength^2)) * (X(a)-Y(b))^2 );
    end
end
%}
%starK = Kstar';
%starK = (1/charlength^2)*(XY).*K;

%{
for a = 1:size(X, 2)
    for b = 1:size(Y, 2)
        starKstar(a,b) = sigmaf*((1/charlength^2) - ((1/charlength^2)*(X(a)-Y(b)))^2 ) * exp( (-1/(2*charlength^2)) * (X(a)-Y(b))^2 );
    end
end
%}

%Cov = ((1/charlength^2)*ones(lengthY,lengthX) - ((1/charlength^2)^2*XY2 )) .* K;
Cov = -(K.*XY2)./(charlength^4) + K./(charlength^4);

end