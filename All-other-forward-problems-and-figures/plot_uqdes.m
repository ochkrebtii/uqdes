function [uensemble,t,logftime] ...
    = plot_uqdes(sspan,N,kernel,lambda,alpha,odefn,odesoln,u0,theta)

%% plot_uqdes.m 
% Run forward probabilistic simulation and plot draws at a selected number
% of iterations
%
% Inputs
% sspan = 1x2 vector containing the upper and lower bound of the domain of
% integration
% N = an integer for the discretizatio mesh size
% odefn = a function handle for the system of ODEs
% u0 = 1xM vector of initial conditions
% theta = parameters used in the system of ODEs, dimension of theta depends
% on the inputs of the function handle odefn
% 
% Outputs
% This function returns 
% ueuler: a NxM vector of approximate ODE solution evaluated at the N grid
% points
% t: the discretization grid at which the function was evaluated (for the
% moment these are equally spaced)
% logftime: gives the log computation time in seconds
% 
% e.g. [umat,tvec,lft] = uqdes([0,10],5,50,'sqexp',1,50,@simpleode,[-1,0],2);


tic % start timer
figure % set up a single figure

M = length(u0); B = 5; 
s = linspace(sspan(1),sspan(2),N); t = s; ds = s(2)-s(1);
if size(u0,1)>size(u0,2); u0 = u0'; end
plotcounter = 0;
% choose colors for five individual draws 
colors = {[0.4, 0.5647058823529412, 0.23529411764705882],...
    [0.6901960784313725,0.3568627450980392,0.792156862745098],...
    [0.7725490196078432,0.3607843137254902,0.19607843137254902],...
    [0.39215686274509803,0.48627450980392156,0.7019607843137254],...
    [0.7647058823529411,0.3058823529411765,0.4823529411764706]};

[ueuler,seuler,logftime] = euler(sspan,N,odefn,u0,theta);

% sort out the kernel to use
switch lower(kernel)
    case {'sqexp'}
        kern = 'se'; %trim = N;
    case {'uniform'}
        kern = 'un'; %trim = 2*ceil(lambda/ds);
    otherwise
        disp('Unknown kernel -- try agaian'); return
end

% define appropriate kernel convolutions
QQ = str2func(strcat('QQ1d_',kern));
RR = str2func(strcat('RR1d_',kern));
QR = str2func(strcat('QR1d_',kern));
uensemble = repmat(u0,[N,1,B]);
u_plot = repmat(u0,[N,1,B]);
f = odefn(s(1),uensemble(1,:,:),theta);
%m_deriv_svec = repmat(f,[N,1,1]);
m_deriv_svec = zeros(N,M,B);
m_state_svec =  uensemble + bsxfun(@times,m_deriv_svec,repmat(s',1,2));
C_deriv_ssmat = RR(s,s,lambda,sspan(1),sspan(2))/alpha;
C_state_ssmat = QQ(s,s,lambda,sspan(1),sspan(2))/alpha;
C_cross1_ssmat = QR(s,s,lambda,sspan(1),sspan(2))/alpha;
kinv = 1/(C_deriv_ssmat(1,1));
f_diff = kinv*(f-m_deriv_svec(1,:,:));
randnNums  = randn(N,M,B);  % generate random numbers outside of loop

%% run algorithm

for n = 0 : N-1
    if n > 0  
        ind = 1:N; endind = 1:N;
        %ind = max(n+1,max(1,n-trim)):min(N,min(N,n+trim));
        %endind = n+1:N;
        m_state_svec(endind,:,:) = m_state_svec(endind,:,:) + bsxfun(@times,C_cross1_ssmat(endind,n),f_diff(1,:,:));
        m_deriv_svec(ind,:,:) = m_deriv_svec(ind,:,:) + bsxfun(@times,C_deriv_ssmat(ind,n),f_diff(1,:,:)); 
        C_state_ssmat(endind,endind) = C_state_ssmat(endind,endind) - kinv*C_cross1_ssmat(endind,n)*(C_cross1_ssmat(endind,n))';
        C_cross1_ssmat(endind,ind) = C_cross1_ssmat(endind,ind) - kinv*C_cross1_ssmat(endind,n)*C_deriv_ssmat(n,ind);
        C_deriv_ssmat(ind,ind) = C_deriv_ssmat(ind,ind) - kinv*C_deriv_ssmat(ind,n)*C_deriv_ssmat(n,ind);
        uensemble(n+1,:,:) = m_state_svec(n+1,:,:) + randnNums(n,:,:)*sqrt(C_state_ssmat(n+1,n+1)); 
        kinv = 1/(C_deriv_ssmat(n+1,n+1)+C_deriv_ssmat(n,n));
        f_diff = kinv*(odefn(s(n+1),uensemble(n+1,:,:),theta) - m_deriv_svec(n+1,:,:));
    end
    
    %% plot subfigure
    if n == 0 || n ==12 || n == 24 || n == 36 || n == 48
        plotcounter = plotcounter+1;
        [Lv,dv] = eig(C_state_ssmat); dv = diag(dv);
        for mm = 1:M
            for bb = 1:B
                u_plot(:,mm,bb) = m_state_svec(:,mm,bb) + real(Lv*(diag(dv.^0.5)))*randnNums(:,mm,bb);
            end
        end
       for mm = 1:M
            f = subaxis(2,5,plotcounter+5*(mm==2));
            hold on
            for bb = 1:B
                plot(s,u_plot(:,mm,bb),'color',colors{bb});
                if n > 0
                    plot(s(n+1),m_state_svec(n+1,mm,bb),'o','markersize',3,'linewidth',1,'color',colors{bb});
                end
                ylims = get(gca,'ylim');
                if (floor(ylims(2))-ceil(ylims(1))) < 6
                    y_label_locations = ceil(ylims(1)):1:floor(ylims(2));
                elseif (floor(ylims(2))-ceil(ylims(1))) < 12
                    if mod(ceil(ylims(1)),2)==0
                        y_label_locations = ceil(ylims(1))+1*(mm==1):2:floor(ylims(2));
                    else
                        y_label_locations = ceil(ylims(1))+1*(mm==2):2:floor(ylims(2));
                    end
                else
                    if mod(ceil(ylims(1)),3)==0
                        y_label_locations = ceil(ylims(1))+1*(mm==1):3:floor(ylims(2));
                    else
                        y_label_locations = (ceil(ylims(1))-mod(ceil(ylims(1)),3)*(mm==2))+1*(mm==1):3:floor(ylims(2));
                    end
                end
                % set(gca,'yTick',y_label_locations)
                if n > 0 && n < N; set(gca,'yTick',[]); end;
            end
            
            exact = odesoln(s,theta);
            plot(s,exact(mm,:),'r--','linewidth',1);
            %plot(s,m_state_svec(1:N,mm,1),'g-','linewidth',1);
            plot(seuler(1:n+1),ueuler(1:n+1,mm),'oc','linewidth',1);
            if mm==1; set(f,'XTickLabel',[]);
            else xlabel(f,'t','Fontsize',18);
            end
            hold off
            box on
            ylims = get(f,'YLim');
            set(f,'Linewidth',1,'Fontsize',18)
            axis([sspan(1),sspan(2),ylims(1),ylims(2)])
       end
    end 
end

logftime = log(toc);    % end timer

end