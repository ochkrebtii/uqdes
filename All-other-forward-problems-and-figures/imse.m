function [IMSE] = imse(t,fn,fntrue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the integrated MSE between fn and fntrue
% using the trapezoidal rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    IMSE = 0;
    
    nsolves = size(fn,3);
    
    for ii = 1:nsolves
       IMSE = IMSE + sum(trapz((fn(:,:,ii)-fntrue).^2));
    end
    
    IMSE = IMSE/nsolves;
    
    if ~isreal(IMSE)
        IMSE = NaN;
    end
    
end