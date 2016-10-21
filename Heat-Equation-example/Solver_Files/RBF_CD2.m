function [ Cov ] = RBF_CD2( hyperparameters, T, DerivT )


sigmaf     = hyperparameters(1);
charlength = hyperparameters(2);

lengthX = size(T, 2);
lengthY = size(DerivT, 2);

%XY  = T(ones(lengthY,1),:)' - DerivT(ones(lengthX,1),:);
XY  = DerivT(ones(lengthX,1),:) - T(ones(lengthY,1),:)';
XY2 = XY.^2;

%{
for a = 1:size(X, 2)
    for b = 1:size(Y, 2)
        K(a,b) = sigmaf*exp( (-1/(2*charlength^2)) * (X(a)-Y(b))^2 );
    end
end
%}
a = (sigmaf./(sqrt(2*pi)*charlength))*exp( (-1/(2*charlength^2)) * XY2 );

%{
for a = 1:size(X, 2)
    for b = 1:size(Y, 2)
        Kstar1(a,b) = -(sigmaf/charlength^2)*(X(a)-Y(b))*exp( (-1/(2*charlength^2)) * (X(a)-Y(b))^2 );
    end
end
%}
b = -(1/charlength^2).*(XY);

Cov = a.*b.*b + a.*(-1/charlength^2);

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

%starKstar = ((1/charlength^2)*ones(lengthY,lengthX) - ((1/charlength^2)^2*XY2 )) .* K;

end