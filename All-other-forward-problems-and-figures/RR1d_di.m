function rrt = RR1d_di(u1,v1,w,a,b)

% u1=[0:0.01:1];
% v1=[0:0.01:1];
% w=0.25;
% a=0;
% b=1;


if (size(u1,1)<size(u1,2))
    u1 = u1';
end

if (size(v1,1)<size(v1,2))
    v1 = v1';
end

u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);

rrt = (w<=(b-a)/2).*(...
      (1./6).*(feval(@(x) x.*(6.*u.*v - 3.*x.*(u+v) + 2.*(x.^2)), min(min(u+w,v+w),b-w))...
           - feval(@(x) x.*(6.*u.*v - 3.*x.*(u+v) + 2.*(x.^2)), max(max(u-w,v-w),a+w))).*(min(min(u+w,v+w),b-w) >= max(max(u-w,v-w),a+w))...
     +(1./4).*(feval(@(x) x.*(a - 2.*u + w).*(a - 2.*v + w) + (x.^2).*(a - u - v + w) + (1./3).*(x.^3) , a + w )...
           - feval(@(x) x.*(a - 2.*u + w).*(a - 2.*v + w) + (x.^2).*(a - u - v + w) + (1./3).*(x.^3) , max(max(u-w,v-w),a))).*( a + w >= max(u-w,v-w))...
     +(1./4).*(feval(@(x) x.*(b - 2.*u - w).*(b - 2.*v - w) + (x.^2).*(b - u - v - w) + (1./3).*(x.^3) , min(min(u+w,v+w),b))...
           - feval(@(x) x.*(b - 2.*u - w).*(b - 2.*v - w) + (x.^2).*(b - u - v - w) + (1./3).*(x.^3) , b - w)).*(min(u+w,v+w) >= b - w)...      
        );     


% mesh(u,v,rrt)


end