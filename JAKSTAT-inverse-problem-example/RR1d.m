function rrt = RR1d(u1,v1,w,a,b)

% a = 0;
% b = 5; 
% N = 100;
% u1 = linspace(a,b,N)
% v1 = u1
% w = 3*(b-a)/N


if (size(u1,1)<size(u1,2))
    u1 = u1';
end

if (size(v1,1)<size(v1,2))
    v1 = v1';
end

u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);


rrt = ((min(u,v)-max(u,v) + 2.*w).*(max(u,v) - min(u,v) < 2.*w));

% mesh(u1,v1,rrt)

end