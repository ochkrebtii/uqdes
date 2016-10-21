%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Diagnostic Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load MCMC data file

clear alpha;  % avoid conflict with transparency

names = {'k_1';'k_2';'k_3';'k_4';'k_5';'k_6';'\tau';'u^{(1)}(0)';'\alpha';'\lambda'};
%obsnames = {'\mathcal{G}^{(1)}(u(t),\theta)';'\mathcal{G}^{(2)}(u(t),\theta)';'\mathcal{G}^{(3)}(u(t),\theta)';'\mathcal{G}^{(4)}(u(t,\theta),\theta)'};
obsnames = {'{G}^{(1)}(u(t),\theta)';'{G}^{(2)}(u(t),\theta)';'{G}^{(3)}(u(t),\theta)';'G^{(4)}(u(t,\theta),\theta)'};
solnames = {'u^{(1)}(t,\theta)';'u^{(2)}(t,\theta)';'u^{(3)}(t,\theta)';'u^{(4)}(t,\theta)'};

pars_true = NaN(size(names));

priors = {@(x) prior_exp(x,priorpars.ratepars(1));...
    @(x) prior_exp(x,priorpars.ratepars(2));...
    @(x) prior_exp(x,priorpars.ratepars(3));...
    @(x) prior_exp(x,priorpars.ratepars(4));...
    @(x) prior_exp(x,priorpars.ratepars(5));...
    @(x) prior_exp(x,priorpars.ratepars(6));...
    @(x) prior_chisq(x,priorpars.taudf);...
    @(x) prior_normal(x,priorpars.hist1mu,priorpars.hist1var);...
    @(x) prior_lognorm(x,priorpars.alphashift(1),priorpars.alphamean(1),priorpars.alphavar(1));...
    @(x) prior_chisq(x,priorpars.lambdasp(1))};...
    


lightgreycol            = [0.6,0.6,0.6];
white                   = [1,1,1];

npars = length(priors);
nstates = 4;
tpn = nchains;
tc = tpn;

lag                             = 20000;
first                           = iter-lag;%floor(iter/2);
last                            = iter-1;
thin                            = 1;

if last > first && first > 1
    span                        = [first:thin:last];
    startkjk                    = floor(first/floor(niter/kjk));
    endkjk                      = floor(last/floor(niter/kjk));
else
    span                        = [1:thin:iter];
    startkjk                    = floor(1/floor(niter/kjk));
    endkjk                      = floor(last/floor(niter/kjk));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1 - traceplots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

for ii = 1:npars
    
    s               = subaxis(3,4,ii);
    cp              = chainpars{tc}(ii,span);
    p1              = plot(1:iter,chainpars{tpn}(ii,1:iter),'b-','Linewidth',1);
    hold on
    p2              = plot(span,chainpars{tpn}(ii,span),'m-','Linewidth',1);
    
    [f,xi]          = ksdensity(cp);
    [fh,xh]         = hist(cp,7);
    xlims           = get(s,'Xlim');
    ylims           = get(s,'Ylim');
    priorgrid       = linspace(ylims(1),ylims(2),100);
        
    plot(xlims(2)*f/max(f),xi,'r--','Linewidth',1)
    plot(xlims(2)*priors{ii}(priorgrid)/max(priors{ii}(priorgrid(~isinf(priors{ii}(priorgrid))))),priorgrid,'k--','Linewidth',1)
    
    if ~isnan(pars_true(ii))
        li = line([0,xlims(2)],[pars_true(ii),pars_true(ii)]);
        set(li,'Color','g','Linestyle','-','Linewidth',1)
    end
    
    axis tight
    axis square
    hold off
    box on
    
    if ii<4
        title(names{ii})
    else
        xlabel(names{ii})
    end
    
end

s = subaxis(3,4,3*4);
plot(1:iter,log_r_bot{tpn}(1:iter),'b-','Linewidth',1)
hold on
plot(span,log_r_bot{tpn}(span),'m-','Linewidth',1)
xlabel('un-normalized log-posterior')
box on

%export_fig 'Traceplots.png'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2 - posterior solution draws
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
solmat                      = cell2mat(chainSol{tpn});
%transmat                    = cell2mat(chainTrans{tpn});
transmat = [];

for kk = startkjk:endkjk
    
    transmat = horzcat(transmat,obstransform(solmat(:,(1+nstates*kk):(nstates+nstates*kk)),chainpars{tpn}(5:6,kk*100+1)));

end

for state = 1:nstates
    
    s = subaxis(2,nstates,nstates+state);
    set(s,'FontSize',14)
    set(s,'linewidth',2,'box','on');%,'xtick',[],'ytick',[],'hittest','off');
    hold on    
    for kk = startkjk:endkjk
        xflip = [dataArray.evalgrid fliplr(dataArray.evalgrid)];
        yflip = [solmat(:,state+nstates*kk)' fliplr(solmat(:,state+nstates*kk)')];
        mins = min(yflip);
        maxs = max(yflip);
        try 
            p = patch(xflip,yflip,'r','EdgeAlpha',1.5/(endkjk-startkjk),'FaceColor','none','Linesmoothing','on','Linewidth',1);
        catch
            p = patch(xflip,yflip,'r','EdgeAlpha',1,'FaceColor','none','Linesmoothing','on','Linewidth',1);
        end    
    end
    xlabel(solnames{state},'fontsize',16);
    hold off
    box on
    axis tight
    axis square
    
    s= subaxis(2,nstates,state);
    set(s,'FontSize',14)
    set(s,'linewidth',2,'box','on');%,'xtick',[],'ytick',[],'hittest','off');
    hold on   
    for kk = 0:length(startkjk:endkjk)-1           
        xflip = [dataArray.evalgrid fliplr(dataArray.evalgrid)];
        yflip = [transmat(:,state+nstates*kk)' fliplr(transmat(:,state+nstates*kk)')];
        mins = min(yflip);
        maxs = max(yflip);
        try 
            p = patch(xflip,yflip,'r','EdgeAlpha',1.5/(endkjk-startkjk),'FaceColor','none','Linesmoothing','on','Linewidth',1);
        catch
            p = patch(xflip,yflip,'r','EdgeAlpha',1,'FaceColor','none','Linesmoothing','on','Linewidth',1);    
        end
    end
    plot(dataArray.time{state},dataArray.obs{state},'rx','Markersize',5,'Linewidth',2)
    hold off
    title(obsnames{state},'fontsize',16);
    box on
    axis tight
    axis square
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3 - correlattion plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set(gcf,'PaperUnits','centimeters')
% %This sets the units of the current figure (gcf = get current figure) on paper to centimeters. 
% xSize = 8; ySize = 8;
% %These are my size variables, width of 8 and a height of 12, will be used a lot later.
% xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% %Additional coordinates to center the figure on A4-paper
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% %This command sets the position and size of the figure on the paper to the desired values.
% set(gcf,'Position',[0 0 xSize*50 ySize*50])
% 

%get(0,'screensize') 

figure%('position',[1 1 1360 768])

xlims =cell(1,npars);

for chainind = 1:npars
    cormat(chainind,:)  = chainpars{tc}(chainind,span);
    cormax(chainind)    = max(cormat(chainind,:));
    cormin(chainind)    = min(cormat(chainind,:));
end

cormax(1) = 3;
cormax(2) = 2.5;
cormax(3) = 0.16;
cormax(4) = 0.25;
cormin(7) = 1;
cormin(8) = 160;
cormax(9) = 0.9e+005;
cormax(10) = 4.5;

timesfactor = ones(1,npars);
timesfactor(3:6) = timesfactor(3:6)/2;
timesfactor(7:8) = timesfactor(7:8)/1.25;
timesfactor(10) = 1.25;

for rr = 1:npars
    for cc = 1:npars
        s = subaxis(npars,npars,(rr-1)*npars+cc,'Spacing',0.001,'SpacingHoriz',0.001,'SpacingVert',0.001,... 
             'Padding',0,'PaddingRight',0,'PaddingLeft',0,'PaddingTop',0,'PaddingBottom',0,... 
             'Margin',0,'MarginRight',0.1,'MarginLeft',0.1,'MarginTop',0.1,'MarginBottom',0.1);         
        set(s,'FontSize',10)
        set(s,'linewidth',1.25,'box','on');%,'xtick',[],'ytick',[],'hittest','off');
        if rr > cc
            %factor          = [(ylims(2)-ylims(1))/5,(xlims(2)-xlims(1))/5];
            %factor          = [1,1];
            densitytemp     = gkde2(cormat([rr,cc],:)');
            densitytemp.h   = densitytemp.h*1.75;
            density         = gkde2(cormat([rr,cc],:)', densitytemp);
            %plot(cormat(cc,:),cormat(rr,:),'Linestyle','none','Marker','o','MarkerEdgeColor',lightgreycol,'MarkerFaceColor',lightgreycol,'Markersize',1);
            %transparentScatter(cormat(cc,:)',cormat(rr,:)',factor,0.2);
            hold on
            [C H] = contour(density.y,density.x,density.pdf);
            set (H, 'LineWidth', 1.25);
            xpad = 0;%(cormax(cc)-cormin(cc))/100;
            ypad = 0;%(cormax(rr)-cormin(rr))/100;
            axis([min(cormin(cc),pars_true(cc)-xpad) max(cormax(cc),pars_true(cc)+xpad) min(cormin(rr),pars_true(rr)-ypad) max(cormax(rr),pars_true(rr)+ypad)]);
            ylims           = get(s,'ylim');
            xlims{rr}       = get(s,'xlim');
            li = line([pars_true(cc),pars_true(cc)],ylims);
            set(li,'Color','g','Linestyle','-','Linewidth',1)
            li = line(xlims{rr},[pars_true(rr),pars_true(rr)]);
            set(li,'Color','g','Linestyle','-','Linewidth',1)
        elseif rr < cc
            set(s,'Visible','off')
        else
            set(s,'XLim',[cormin(rr) cormax(rr)]);
            %set(s,'XLim',xlims{rr});
            cp              = chainpars{tc}(rr,span);
            [f,xi]          = ksdensity(cp);
            %[fh,xh]         = hist(cp,7);
            %bar(xh,fh/max(fh),'FaceColor',lightgreycol,'EdgeColor',lightgreycol);
            plot(xi,f/max(f),'Linewidth',1,'Color','r','Linestyle','--')
            hold on
            plot(xi,timesfactor(rr)*priors{rr}(xi)/max(priors{rr}(xi)),'Linewidth',1,'Color','k','Linestyle','--')
            ylims = [0,1];
            li = line([pars_true(rr),pars_true(rr)],ylims);
            set(li,'Color','g','Linestyle','-','Linewidth',1)
            xpad = (cormax(cc)-cormin(cc))/20;
            axis([min(cormin(cc),pars_true(cc)-xpad) max(cormax(cc),pars_true(cc)+xpad) ylims(1) ylims(2)+0.25]);
        end
        
        
        
        if rr == npars
            xlabel(names{cc},'fontsize',16);
        end
        
        if cc == 1
            ylabel(names{rr},'fontsize',16);
        end
        
        if cc > 1
            set(s, 'YTickLabel', [], 'YTick', []);
        
        end
        
        set(s,'box','on');

    end
end

%export_fig 'jakstat-corrplots.pdf'

