function [ CovMatrix ] = Covariance_PartialDerivatives_RBF( hyperparameters, X, Y )


sigmaf     = hyperparameters(1);
charlength = hyperparameters(2);

lengthX = size(X, 2);
lengthY = size(Y, 2);

XY2 = ( X(ones(lengthY,1),:)' - Y(ones(lengthX,1),:)).^2;

% PD wrt hyperparameter 1
%{
for a = 1:size(X, 2)
    for b = 1:size(Y, 2)
        CovMatrix{1}(a,b) = exp( (-1/(2*charlength^2)) * (X(a)-Y(b))^2 );
    end
end
%}
CovMatrix{1} = (1/(sqrt(2*pi)*charlength))*exp( (-1/(2*charlength^2)) * XY2 );

% PD wrt hyperparameter 2
%{
for a = 1:size(X, 2)
    for b = 1:size(Y, 2)
        CovMatrix{2}(a,b) = sigmaf * exp( (-1/(2*charlength^2))*(X(a)-Y(b))^2 ) * ( (1/(charlength^3)) * (X(a)-Y(b))^2 );
    end
end
%}
%CovMatrix{2} = sigmaf * exp( (-1/(2*charlength^2))*XY2 ) .* ( (1/(charlength^3)) * XY2 );
CovMatrix{2} = -(sigmaf/(sqrt(2*pi)*charlength^2))*exp( (-1/(2*charlength^2)) * XY2 ) + (sigmaf.*CovMatrix{1}) .* ( (1/(charlength^3)) * XY2 );

% PD wrt hyperparameter 3
CovMatrix{3} = eye(lengthX);


end