function newdata = convertYTcell(datavec,lag)

%% INPUT
% lag                   - length of each lag interval
% datavec               - structure array ordered by observation process (datavec.time, datavec.obs)

%% RETURNS
% newdata               - structure array ordered by time lag (newdata.time, newdata.obs, newdata.n)

a                               = datavec.tspan(2);

nlags                           = floor(range(datavec.tspan(2:3))/lag);

if  floor(range(datavec.tspan(2:3))/lag) == range(datavec.tspan(2:3))/lag
    nlags                       = nlags-1;
end

newdata                         = datavec;
newdata.nlags                   = nlags;

end


