%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Diagnostic Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {'\theta_1';'\theta_2';'\theta_3';'u_0';'v_0';'\tau_1^2';'\tau_2^2';'\alpha';'\beta'; '\lambda_1';'\lambda_2';'L_1';'L_2'};

% SET SPAN
lag=1000;
try
    first=restart_at_iter;
catch
    first=1;
end
first = 10000;
last=iter;
thin=1;
if(last>first) span=[first:thin:last];
else span=[iter-lag:iter];end


% CORRELATION FIGURES
% model parameters
counter = 1;
figure('Color',[1 1 1])
subplot('Position',[0.1 0.1 0.5 0.5]);
for i=1:8
    for j=1:8
        if(i>j)
            subaxis(8,8,(i-1)*8 + j)
            plot(chainpars{tc}(j,span),chainpars{tc}(i,span),...
                'LineStyle','none','Marker','o','MarkerEdgeColor',[0.2 0.2 0.2],...
                'Markersize',1.5)
            if (j==1)
                ylabel(names(i))
            end
            if(i==5)
                xlabel(names(j))
            end
        elseif(i==j)
            subaxis(8,8,(i-1)*8 + j)
            cp = chainpars{tc}(i,span);
            [f,xi] = ksdensity(cp);
            [fh,xh] = hist(cp,7);
            hold on;
            %camroll(270);
            bar(xh,fh/(sum(fh)),'FaceColor',[0.4,0.4,0.4],'EdgeColor',[0.3,0.3,0.3]);
            title(names(i))
            if(i==1)
                ylabel(names(i))
            end
            counter=counter+1;
            box on
            hold off;
        end
    end
end
print -depsc -noui modparscorrs.eps


%auxiliary parameters
figure('Color',[1 1 1])
for i=8:11
    for j=8:11
        if(i>j)
            subaxis(4,4,(i-7-1)*4 + j-7)
            plot(chainpars{tc}(j,span),chainpars{tc}(i,span),...
                'LineStyle','none','Marker','o','MarkerEdgeColor',[0.2 0.2 0.2],...
                'Markersize',1.5)
            if (j==1)
                ylabel(names(i))
            end
            if(i==5)
                xlabel(names(j))
            end
        elseif(i==j)
            subaxis(4,4,(i-7-1)*4 + j-7)
            cp = chainpars{tc}(i,span);
            [f,xi] = ksdensity(cp);
            [fh,xh] = hist(cp,7);
            hold on;
            %camroll(270);
            bar(xh,fh/(sum(fh)),'FaceColor',[0.4,0.4,0.4],'EdgeColor',[0.3,0.3,0.3]);
            title(names(i))
            if(i==1)
                ylabel(names(i))
            end
            counter=counter+1;
            box on
            hold off;
        end
    end
end
print -depsc -noui auxparscorrs.eps



load('MCMC solve the DE each time_run on date_31-Oct-2011.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Diagnostic Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = {'\theta_1';'\theta_2';'\theta_3';'u_0';'v_0';'\tau_1^2';'\tau_2^2'};

% SET SPAN
lag=1000;
try
    first=restart_at_iter;
catch
    first=1;
end
last=iter;
thin=1;
if(last>first) span=[first:thin:last];
else span=[iter-lag:iter];end


% CORRELATION FIGURES
% model parameters
counter = 1;
figure('Color',[1 1 1])
subplot('Position',[0.1 0.1 0.5 0.5]);
for i=1:5
    for j=1:5
        if(i>j)
            subaxis(5,5,(i-1)*5 + j)
            plot(chainpars(j,span),chainpars(i,span),...
                'LineStyle','none','Marker','o','MarkerEdgeColor',[0.2 0.2 0.2],...
                'Markersize',1.5)
            if (j==1)
                ylabel(names(i))
            end
            if(i==5)
                xlabel(names(j))
            end
        elseif(i==j)
            subaxis(5,5,(i-1)*5 + j)
            cp = chainpars(i,span);
            [f,xi] = ksdensity(cp);
            [fh,xh] = hist(cp,7);
            hold on;
            bar(xh,fh/(sum(fh)),'FaceColor',[0.4,0.4,0.4],'EdgeColor',[0.3,0.3,0.3]);
            title(names(i))
            if(i==1)
                ylabel(names(i))
            end
            counter=counter+1;
            box on
            hold off;
        end
    end
end
print -depsc -noui modparscorrs_num.eps

