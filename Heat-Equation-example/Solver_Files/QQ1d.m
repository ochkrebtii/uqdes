function qqt = QQ1d(u1,v1,w,a,b)


% u1=[0:0.1:1]
% v1=[0:0.1:1]
% w1=0.2
% w2=0.1
% a=0
% b=1

if (size(u1,1)<size(u1,2))
    u1 = u1';
end

if (size(v1,1)<size(v1,2))
    v1 = v1';
end

u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);


qqt=((4.*w.^2).*(min(u-w,v-w)-(a+w)).*(min(u-w,v-w)>(a+w))...
    +(2.*w).*((v+w).*min(u-w,v+w) - 0.5.*min(u-w,v+w).^2-(v+w).*max(a+w,v-w) + 0.5.*max(a+w,v-w).^2).*(min(u-w,v+w)>max(a+w,v-w))...
    +((1/3).*min(a+w,min(u-w,v-w)).^3 + (w-a).*min(a+w,min(u-w,v-w)).^2 + (w-a).^2.*min(a+w,min(u-w,v-w))-(1/3).*(a-w).^3 - (w-a).*(a-w).^2 - (w-a).^2.*(a-w)).*(min(a+w,min(u-w,v-w))>(a-w))...
    +(v-a).*(0.5.*min(a+w,u-w).^2 + (w-a).*min(a+w,u-w)- 0.5.*max(a-w,v-w).^2 - (w-a).*max(a-w,v-w)).*(min(a+w,u-w)>max(a-w,v-w))...
    +(2.*w).*((u+w).*min(u+w,v-w)-0.5.*min(u+w,v-w).^2-(u+w).*max(a+w,u-w) + 0.5.*max(a+w,u-w).^2).*(min(u+w,v-w)>max(a+w,u-w))...
    +((u+w).*(v+w).*min(u+w,v+w) - 0.5.*(u+v+2.*w).*min(u+w,v+w).^2 + (1/3).*min(u+w,v+w).^3 - (u+w).*(v+w).*max(a+w,max(u-w,v-w)) + 0.5.*(u+v+2.*w).*max(a+w,max(u-w,v-w)).^2 - (1/3).*max(a+w,max(u-w,v-w)).^3).*(min(u+w,v+w)>max(a+w,max(u-w,v-w)))...
    +(u-a).*(0.5.*min(a+w,v-w).^2 + (w-a).*min(a+w,v-w)- 0.5.*max(a-w,u-w).^2 - (w-a).*max(a-w,u-w)).*(min(a+w,v-w)>max(a-w,u-w))...
    +(u-a).*(v-a).*((a+w)-max(u-w,v-w)).*((a+w)>max(u-w,v-w))...
    );
        


end
