function [sampled_function] = Example_Lane_Emden_MBVP(pars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example_Lane_Emden_MBVP.m
% This function samples realizations from the forward probabilistic
% solution of a mixed boundary value problem as described in the paper 
% Bayesian Solution Uncertainty Quantification for Differential Equations
% by O.A. Chkrebtii, D.A. Campbell, B. Calderhead, M.A. Girolami.

% Inputs
% pars = model parameter(s) 
% 
% Outputs
% this program outputs a mat file containing a number (nchains) of parallel
% chains, within the variable chainpars.  The last chain is obtained with a
% likelihood tempered by one.
% 
% e.g. sampled_function = Example_Lane_Emden_MBVP([1,2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('This will take some time depending on the number of MCMC iterations and the number of compute cores used')

% number of MCMC iterations
niter       = 2000; 
% number of proposal variance updates (note: stop updating half way through)
kjk         = niter/50;  
% number of solver knots
Nknots      = 100;
% number of parallel chains
nchains      = 6;


% set random seed (or comment out to generate one automatically)
%RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))


%%%%%%%%%%%%%%%%%%%%%%%
% Set true functions
%%%%%%%%%%%%%%%%%%%%%%%

% Lorenz system
odefn = @emdenmbv;
LogPostProb = @problpp;


%%%%%%%%%%%%%%%%%%%%%%%
% Set solver specifications
%%%%%%%%%%%%%%%%%%%%%%%

% domain of integration
a                   = 0.5;
b                   = 1;
tspan               = [a b];

% number of solver knots
if isempty(Nknots)
    Nknots  = 100;
end


% solver grid points
tgrid      = linspace(tspan(1),tspan(2),Nknots);


%%%%%%%%%%%%%%%%%%%%%%%
% Set true parameters
%%%%%%%%%%%%%%%%%%%%%%%

% model parameters
if nargin==0 || isempty(pars)
    pars                = [1,2];
end
xfirst              = [ NaN ; -(tspan(1)/3).*(1 + (tspan(1).^2)/3).^(-3/2) ];
xlast               = [ (1 + (tspan(2).^2)/3).^(-1/2) ; NaN ];

alpha               = 1;
lambda              = 2*range(tspan)/(Nknots-1);

% data
data = [xlast(1),NaN];
time = tspan(2);

if sum(tgrid==tspan(2))==0
    tspan = [tgrid,tspan(2)];
end

% error term variances
%sigma_true      = [1,1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the MCMC stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  name the output file
filename    = strcat('Emden_prob_soln');
filename    = strcat(filename,'_',date,'.mat');     % mcmc runs are saved in this file eventually


%% set up temperatures for parallel tempering ( set tp = 1 for regular MCMC )
tp              = linspace(0.1,1,nchains);

%% Prior Hyperparameters

% hyperparameters for model parameters
mu_start        = 1.5;
var_start       = (2*(xlast(1) - xfirst(2)))^2;


% put all these together

prior_pars  =   [mu_start,var_start];

%% Starting values

% load starting values from this file, or define them here
new_stepvar_chain = 0.05;

% empty starting value for solution realization
x_start = [];
pars_start = [pars,[mu_start,xfirst(2)],xlast(1),alpha,lambda]; %,sigma_true];


% evaluate and display log posteriors for first iteration
disp('---------------------')
disp('initial un-normalized log posteriors')
for jj=1:length(tp)
    stepvar_chain{jj}=new_stepvar_chain/tp(jj);
    [log_r_bot1(jj),X{jj}] = LogPostProb(Nknots,tspan,tgrid,pars_start,x_start,time,data,prior_pars,odefn,tp(jj));
    disp(log_r_bot1(jj))
end

disp('---------------------')

for j=1:kjk
    swap_counter{j} = zeros(1,length(tp));
end

%% Chain stuff

kjk_stop_fix = kjk;
k=1;
div_by = floor(niter/kjk);
% log denominator matrix

for jj=1:length(tp)
    rate{jj}             = zeros(kjk+1,1);                                          % keep track of the acceptance rate
    accepts{jj}          = zeros(kjk+1,1);                                         % count when a proposal is accepted
    log_r_bot{jj}        = [log_r_bot1(jj),zeros(1,niter)];
    chainpars{jj}        = [pars_start',zeros(length(pars_start),niter)];           % parameter matrix
    chainxs{jj}          = [X{jj},zeros(length(X{jj}),2*ceil(niter/(div_by)))];
end;

iter = 1;                                                % start iteration counter

itertime=0;

% save filename
save(filename)



%% FINALLY RUN THE MCMC CHAIN

for cc=1:length(tp)
    tempchainpars{cc}=[];
    templog_r_bot{cc}=[];
    tempaccepts{cc}=[];
    temptime{cc}=time;
    temptgrid{cc}=tgrid;
    temptspan{cc}=tspan;
    tempdata{cc}=data;
end

for iter = 1:niter
    tic
    rrr = iter+1;
    
    % Propose to swap between chains.
    chain1 = ceil(rand*(length(tp)+2));
    chain2 = ceil(rand*(length(tp)+2));
    chains_swap=sort([chain1,chain2]);
    if(chains_swap(1)==chains_swap(2) ||chains_swap(1)<1||chains_swap(2)>length(chainpars))
        chains_swap=[0,0];
    else
        prop_pars{1}=chainpars{chains_swap(2)}(:,iter);
        prop_pars{2}=chainpars{chains_swap(1)}(:,iter);
        prop_X{1}=X{chains_swap(2)};
        prop_X{2}=X{chains_swap(1)};
    end
    rr(iter,1:2)=chains_swap;
    
    % if there is a swap use these:
    if chains_swap(1)>0
        
        [t_newpars,t_newX,lrb,s_count] = swap_pars_emden(Nknots,tspan,tgrid,...
            {chainpars{chains_swap(1)}(:,iter),chainpars{chains_swap(2)}(:,iter)},...
            {X{chains_swap(1)},X{chains_swap(2)}},...
            [log_r_bot{chains_swap(1)}(iter),log_r_bot{chains_swap(2)}(iter)],...
            prop_pars,prop_X,prior_pars,time,data,odefn,0,LogPostProb,tp(chains_swap));
        
        log_r_bot{chains_swap(1)}(rrr) = lrb(1);
        log_r_bot{chains_swap(2)}(rrr) = lrb(2);
        X{chains_swap(1)} = t_newX{1};
        X{chains_swap(2)} = t_newX{2};
        chainpars{chains_swap(1)}(:,rrr) = t_newpars{1};
        chainpars{chains_swap(2)}(:,rrr) = t_newpars{2};
        
        swap_counter{k}([chains_swap(1),chains_swap(2)])...
            = swap_counter{k}([chains_swap(1),chains_swap(2)]) + s_count;
        rr(iter,3)=s_count;
        
        
        
        % if a swap is proposed
        parfor(tc=1:length(tp)) %non-swap chains only
            %temptime=time;
            if(tc~=chains_swap(1)&&tc~=chains_swap(2)&& tc<length(tp)+1)
                [tempchainpars{tc},X{tc}, templog_r_bot{tc}, tempaccepts{tc},] =...
                    sampler(Nknots,temptspan{tc},temptgrid{tc},chainpars{tc}(:,iter),X{tc},log_r_bot{tc}(iter),...
                    prior_pars,temptime{tc},tempdata{tc},accepts{tc}(k,:),stepvar_chain{tc},odefn,LogPostProb,tp(tc));
            end
        end
        
        for(tc=1:length(tp))
            if(tc~=chains_swap(1)&&tc~=chains_swap(2)&& tc<length(tp)+1)
                chainpars{tc}(:,rrr) = tempchainpars{tc};
                log_r_bot{tc}(rrr) = templog_r_bot{tc};
                accepts{tc}(k,:) = tempaccepts{tc};
            end
        end
        
    else    % If no swap is proposed, then use these:
        
        parfor(tc=1:length(tp)) % non-swap chains only
            if(tc~=chains_swap(1)&&tc~=chains_swap(2)&& tc<length(tp)+1)
                [tempchainpars{tc},X{tc},templog_r_bot{tc},tempaccepts{tc}] =...
                    sampler(Nknots,temptspan{tc},temptgrid{tc},chainpars{tc}(:,iter),X{tc},log_r_bot{tc}(iter),...
                    prior_pars,temptime{tc},tempdata{tc},accepts{tc}(k,:),stepvar_chain{tc},odefn,LogPostProb,tp(tc));
            end
        end
        
        for(tc=1:length(tp))
            %if(tc~=chains_swap(1)&&tc~=chains_swap(2)&& tc<length(tp)+1)
            if(tc<length(tp)+1)
                chainpars{tc}(:,rrr) = tempchainpars{tc};
                log_r_bot{tc}(rrr) = templog_r_bot{tc};
                accepts{tc}(k,:) = tempaccepts{tc};
            end
        end
        
    end
    
    
    if rem(iter,div_by)==0
        
        disp(strcat('number of iterations so far: ',num2str(iter)))
        disp(strcat('minutes since start: ',num2str(itertime/60)))
        disp('-----------------')
        %if(rem(iter,div_by*10)==0)  % record sampled stuff every 5000 iterations
        for(mm=1:length(tp))
            chainxs{mm}(:,2*iter/(div_by)+1:2*iter/(div_by)+2) = X{mm};
            disp(accepts{mm}(k))
        end
        %end
        
        if iter<=niter/2
            
            disp('-----------------')
            
            for mm=1:length(tp)
                rate{mm}(k,:)=(1+accepts{mm}(k,:))./...         % add one accept and one non-accept before computing the ratio.
                    (div_by-sum(sum(rr((k-1)*div_by+1:k*div_by,1:2)==mm))+2);
                %disp(rate{mm}(k,:))
                
                if ((rate{mm}(k,1)<.34 || rate{mm}(k,1)>.44)&& k<=kjk_stop_fix)  % Note that I stop making adjustments by kjk_stop_fix iterations
                    stepvar_chain{mm}=stepvar_chain{mm}*rate{mm}(k,1)/.39;
                end
                % target acceptance rates are 44% for one dimension decaying down to 23% for 5 and more dimensions
            end
            
        end
        
        save(filename)
        
        k=k+1;
        
        
    end
    
    itertime=itertime+toc;
end

disp(strcat('total time elapsed: ',num2str(itertime/60)))

% return sampled function
tempind = randi([floor(iter/(2*div_by)),iter/div_by]);
sampled_function = chainxs{length(tp)}(:,2*tempind+1:2*tempind+2);

save(filename)
make_diagnostic_plots

end
