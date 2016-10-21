function rrt = RR1d_se(u1,v1,w,a,b)

% u1=[0:0.1:1]
% v1=[0:0.1:1]
% w=0.2
% a{1} = zeros(size(u1));
% a{2} = zeros(size(v1));

if (size(u1,1)<size(u1,2))
    u1 = u1';
end

if (size(v1,1)<size(v1,2))
    v1 = v1';
end

u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);

%rrt = exp(-((u-v).^2)/(4.*w.^2))/(2.*sqrt(pi).*w);
rrt = exp(-((u-v).^2)/(4.*w.^2))*sqrt(pi).*w;

end