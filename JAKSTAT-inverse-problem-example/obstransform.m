function [obsvec] ...
    = obstransform(statevec,scalingfactors)

obsvec = statevec;

obsvec(:,1)     = scalingfactors(1)*(statevec(:,2) + 2*statevec(:,3));
obsvec(:,2)     = scalingfactors(2)*(statevec(:,1) + statevec(:,2) + 2*statevec(:,3));
obsvec(:,3)     = statevec(:,1);
obsvec(:,4)     = statevec(:,3)./(statevec(:,2)+statevec(:,3));

   
end
