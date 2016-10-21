function r = emdenmbv(t,y,pars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% emdenmbv
% p is a 2x1 vector of parameters, p=(1,2) if empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



try
    
    if (isempty(pars))
        pars=[1;2];
    end
    r               =     y ;
    r(1)            =     pars(1)*y(2) ;
    r(2)            =  - (pars(2)/t).*y(2) - y(1).^5 ;
    
catch
    
    r = [];
    
end


end