function rrt = RR1d_un(u1,v1,lambda,a,b)

if (size(u1,1)<size(u1,2))
    u1 = u1';
end

if (size(v1,1)<size(v1,2))
    v1 = v1';
end

u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);


rrt = ((min(u,v)-max(u,v) + 2.*lambda).*(min(u,v) - max(u,v) > -2.*lambda));

% mesh(u1,v1,rrt)

end