function [past sv sd] = jakstathistory(t,pars)

if isempty(pars) 
    past = zeros(1,4); 
else
    past = pars';
end

sv = [0;0;0;0]';
sd = [0;0;0;0]';

end