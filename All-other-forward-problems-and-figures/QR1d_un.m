function qrt = QR1d_un(u1,v1,lambda,a,b)

if (size(u1,1)<size(u1,2))
    u1 = u1';
end

if (size(v1,1)<size(v1,2))
    v1 = v1';
end

u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);

qrt=(...
     (2.*lambda).*(min(u-lambda,v+lambda) - max(a+lambda,v-lambda)).*(min(u-lambda,v+lambda) > max(a+lambda,v-lambda))...
    +(0.5.*min(a+lambda,min(u-lambda,v+lambda)).^2 + (lambda-a).*min(a+lambda,min(u-lambda,v+lambda)) - 0.5.*(v-lambda).^2 -(lambda-a).*(v-lambda)).*(min(a+lambda,min(u-lambda,v+lambda)) > (v-lambda))...
    +((u+lambda).*min(u+lambda,v+lambda) - 0.5.*min(u+lambda,v+lambda).^2 - (u+lambda).*max(a+lambda,max(u-lambda,v-lambda)) + 0.5.*max(a+lambda,max(u-lambda,v-lambda)).^2).*(min(u+lambda,v+lambda) > max(a+lambda,max(u-lambda,v-lambda)))...
    + (u-a).*(-max(u,v) + a + 2.*lambda).*(a-max(u,v) > -2.*lambda) ...
    );


end