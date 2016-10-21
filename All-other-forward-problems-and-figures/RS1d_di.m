function rst = RS1d_di(u1,v1,w,a,b)

% u1=[0:0.01:0.3];
% v1=[0:0.01:0.3];
% w=0.1;
% a=0;
% b=0.3;


if (size(u1,1)<size(u1,2))
    u1 = u1';
end

if (size(v1,1)<size(v1,2))
    v1 = v1';
end

u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);

rst = (w<=(b-a)/2).*(...
      -(feval(@(x) (1./2).*(x.^2) - u.*x, min(min(u+w,v+w),b-w))...
           - feval(@(x) (1./2).*(x.^2) - u.*x, max(max(u-w,v-w),a+w))).*(min(min(u+w,v+w),b-w) >= max(max(u-w,v-w),a+w))...
     -(1./2).*(feval(@(x) x.*(a - 2.*u + w) + (1./2).*(x.^2) , min(a+w,b-w))...
           - feval(@(x) x.*(a - 2.*u + w) + (1./2).*(x.^2) , max(max(u-w,v-w),a))).*(min(a+w,b-w) >= max(max(u-w,v-w),a))...
     -(1./2).*(feval(@(x) x.*(b - 2.*u - w) + (1./2).*(x.^2) , min(min(u+w,v+w),b))...
           - feval(@(x) x.*(b - 2.*u - w) + (1./2).*(x.^2) , max(a+w,b-w))).*(min(min(u+w,v+w),b) >= max(a+w,b-w))...      
     );     


% mesh(u,v,rst)

end