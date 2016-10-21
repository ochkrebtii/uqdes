function [ Cov ] = RBF_CC( hyperparameters, T1, T2 )


sigmaf     = hyperparameters(1);
charlength = hyperparameters(2);

lengthT1 = size(T1, 2);
lengthT2 = size(T2, 2);

%{
CovMatrix = zeros(lengthT1,lengthT2);
for a = 1:lengthT1
    for b = 1:lengthT2
        CovMatrix(a,b) = sigmaf*exp( (-1/(2*charlength^2)) * (T1(a)-T2(b))^2 );
    end
end
%}

% Vectorized
%CovMatrix = sigmaf*exp( (-1/(2*charlength^2))*(repmat(X',1,lengthY) - repmat(Y,lengthX,1)).^2 );

% Faster than repmat
T1T2Sq = ( T2(ones(lengthT1,1),:) - T1(ones(lengthT2,1),:)' ).^2;
Cov    = (sigmaf/(sqrt(2*pi)*charlength))*exp( (-1/(2*charlength^2))* T1T2Sq );


end