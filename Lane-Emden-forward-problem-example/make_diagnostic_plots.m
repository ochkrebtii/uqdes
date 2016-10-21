%% Make Diagnostic Plots

par = 3;
tempind = floor(iter/(2*div_by));

figure
subplot(1,3,1)
hold on
for kk = 2*tempind+1:2:2*iter/(div_by)+1
    xrealization{kk}    = chainxs{nchains}(:,kk);
    xflip               = [tgrid fliplr(tgrid)];
    yflip               = [xrealization{kk}' fliplr(xrealization{kk}')];
    p                   = patch(xflip,yflip,'r','EdgeAlpha',0.2,'FaceColor','none');
end
plot(time,data(:,1),'.')
box on
hold off
axis([0.5,1.001,0.8,2])
ylabel('u')
xlabel('t')
set(gca,'Linewidth',1.5,'fontsize',26);

subplot(1,3,2)
hold on
for kk = 2*tempind+2:2:2*iter/(div_by)+2
    xrealization{kk}    = chainxs{nchains}(:,kk);
    xflip               = [tgrid fliplr(tgrid)];
    yflip               = [xrealization{kk}' fliplr(xrealization{kk}')];
    p                   = patch(xflip,yflip,'r','EdgeAlpha',0.2,'FaceColor','none');
end
plot(time,data(:,2),'.')
box on
hold off
axis([0.5,1.001,-3,0])
ylabel('v')
xlabel('t')
set(gca,'Linewidth',1.5,'fontsize',26);



%% histograms

subplot(1,3,3)

cp = chainpars{tc}(3,1000:2:end);
cp1 = cp(cp<=1.05 & cp>=0.925);
cp2 = cp(cp<=1.95 & cp>=1.85);
[f1,xi1] = ksdensity(cp1);
[fh1,xh1] = hist(cp1,12);
[f2,xi2] = ksdensity(cp2);
[fh2,xh2] = hist(cp2,13);

f = [f1/sum(f1),f2/sum(f2)];
y = 6*f;
x = [xi1,xi2];

fh = [fh1/sum(fh1),fh2/sum(fh2)];
z = fh;
w = [xh1,xh2];

% define plotting parameters
start = 1.02;
stop = 1.84;
width = 0;


% erase unused data
y(x>start & x<stop)=[];
x(x>start & x<stop)=[];

z(w>start & w<stop)=[];
w(w>start & w<stop)=[];

% map to new xaxis, leaving a space 'width' wide
x2=x;
x2(x2>=stop)=x2(x2>=stop)-(stop-start-width);

w2=w;
w2(w2>=stop)=w2(w2>=stop)-(stop-start-width);

%h2=bar(w2,z,'FaceColor',[0.4,0.4,0.4],'EdgeColor',[0.4,0.4,0.4]);
hold on
h1=plot(x2,y,'b-','Linewidth',1);
hold off

axis([min(x2)-0.01,max(x2)+0.01,0,0.2])
box on

ytick=get(gca,'YTick');
t1=text(start+width/2,ytick(1),'//','fontsize',26);
t2=text(start+width/2,ytick(max(length(ytick))),'//','fontsize',26);
% For y-axis breaks, use set(t1,'rotation',270);

% remap tick marks, and 'erase' them in the gap
xtick=get(gca,'XTick');
dtick=xtick(2)-xtick(1);
gap=floor(width/dtick);
last=max(xtick(xtick<=start)); % last tick mark in LH dataset
next=min(xtick(xtick>=(last+dtick*(1+gap)))); % first tick mark within RH dataset
offset=size(x2(x2>last&x2<next),2)*(x(2)-x(1));

for i=1:sum(xtick>(last+gap))
    xtick(find(xtick==last)+i+gap)=stop+offset+dtick*(i-1);
end
    
for i=1:length(xtick)
    if xtick(i)>last & xtick(i)<next
        xticklabel{i}=sprintf('%.2f',[]);
    else
        xticklabel{i}=sprintf('%.2f',xtick(i));
    end
end;

xlabel('u(t=0.5)')
box on
ylabel('un-normalized log likelihood')
set(gca,'xticklabel',xticklabel,'Linewidth',1.5);
set(gca,'yticklabel',[]);
set(gca,'Linewidth',1.5,'fontsize',26);




