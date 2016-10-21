function [] = JAKSTAT_inference(niter,kjk,N,filename,newcovsfilename,restart_at_iter)

%% JAKSTAT_inference
% This function incorporates a probabilistic delay differentail equation
% solver into a parallel tempering MCMC algorithm for inference on the
% JAK-STAT system of protein network dynamics.
% This method is described in the paper Bayesian Uncertainty Quantification
% for Differential Equations, by O.A. Chkrebtii, D.A. Campbell, M.A.
% Girolami, B. Calderhead.

% INPUTS
% niter = number of iterations for the MCMC
% kjk = total number of displayed summaries produced (covariance updates stop after 3*niter/kjk iterations) 
% N = size of discretization grid for the solver
% filename = name of the output file
% newcovsfilename = filename of tuned transition covariance matrix
% restart_at_iter = iteration at which to restart MCMC after it was stopped
% (make sure restart_at_iter<=iter) and that the file filename.mat is in the work ing folder 
% 
% OUTPUTS
% this program outputs a mat file containing a number (nchains) of parallel
% chains, within the variable chainpars.  The last chain is obtained with a
% likelihood tempered by one.
%
% e.g. JAKSTAT_inference(1000000,10000,500) runs an MCMC with 1000000 iterations and a discretization grid size of 500 (recommendation: use N >= 500 for this system)
% We recommend using multiple parallel cores for this function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET UP THE MCMC

if nargin == 3
   
    %% set up nchains temperatures for parallel tempering ( set tp = 1 for single chain )
    nchains = 10;
    tp = linspace(0.5,1,nchains);
    
    disp(['You are using ',num2str(nchains),' parallel chain(s)'])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Set random seed
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %RandStream.setGlobalStream(RandStream('mt19937ar','seed',83))
    
    %disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    %disp([rand,normrnd(0,1,1,1)])
    %disp(' remember that you are using a fixed seed here')
    %disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Set true functions
    %%%%%%%%%%%%%%%%%%%%%%%
    
    ddefn = @jakstatfundde;       % RHS function
    ddehist = @jakstathistory;      % initial function
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Set starting model & auxiliary parameters
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %load EPO data
    load('JakStatEpo')
    
    %model                   = [1.7294,2.5297,0.098385,0.19643,0.0064224,0.0046532]';  % parameter values
    %histpars                = [207.1,0,0,0]';
    
    dataArray.EpoRAobs      = pEpoR;
    dataArray.EpoRAtime     = time;
    
    %laglength               = 4;                      % lag length
    
    a                       = 0;                        % interval lower bound
    b                       = 60;                       % interval upper bound
    
    %lambda                  = 3*(b-a)/(N-1);            % length-scale
    %alpha                   = 1000;                    % prior precision should increase with the grid size
    %epoauxpars              = [lambda,Inf]';
    
    %pars_start              = vertcat(model,laglength,histpars(1),alpha,lambda);
    %pars_start              = vertcat(model,laglength,histpars(1),alpha,lambda,epoauxpars);
    
    % load starting values from this file, or define them here
    
    %new_stepvar_chain       = diag([ones(1,7)*0.0001,0.01,100,0.01]);
    
    %load('new_pilot_cov_pars');
    
    %pars_start              = vertcat(model,laglength,histpars(1),alpha,lambda);
    
    %if nchains == 10
    %    load('loadcovpars.mat');
    %else
    
    pars_start = [ 1.9535,0.90475,0.089,0.17362,0.0060814,0.0044625,3.2553,215.25,15581,2.1649]';
    new_stepvar_chain = diag([0.01,0.01,0.0001,0.0001,1e-8,1e-8,0.01,20,1e8,0.1]);
    
    %end
    
    %PropSoln_start          = [];                       % empty starting value for solution realization
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load simulated data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load('JakStatData');
    
    dataArray.time = Tcell;
    dataArray.obs = Ycell;
    dataArray.sigsq = Sigsqcell;
    dataArray.n = cellfun(@(u) length(u),Tcell,'Uni',0);
    dataArray.stateheaders = stateheaders;      % array with model state names
    dataArray.obsheaders = obsheaders;        % array with observation state names
    dataArray.nstates = length(dataArray.n);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Set other solver specifications
    %%%%%%%%%%%%%%%%%%%%%%%
    
    N                       = 500;                   % number of grid knots
    npars                   = length(pars_start);
    tspan                   = [a b];                 % interval range
    
    nstates                 = length(ddehist(a,[]));
    nevalpts                = N;                      % number of grid points for evaluation and plotting
    evalgrid                = unique(sort(horzcat(linspace(a,b,nevalpts),cell2mat(Tcell))));
    
    dataArray.evalgrid      = evalgrid;
    dataArray.tspan         = tspan;
    
    % obtain indices of data locations within evalgrid
    for state = 1:nstates
        
        [~,dataArray.inds{state}] = intersect(evalgrid,dataArray.time{state});
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %% Prior Hyperparameters
    %%%%%%%%%%%%%%%%%%%%%%%
    
    priorpars.ratepars = [1,1,1,1,1,1];
    priorpars.taudf = 6;
    priorpars.hist1mu = Ycell{3};  %EB
    priorpars.hist1var = 40^2;
    priorpars.alphashift = 100;
    priorpars.alphamean = 10;
    priorpars.alphavar = 1;
    priorpars.lambdasp = 1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the MCMC stuff
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    temp = clock;
    
    filename = strcat('output_DDE_ParEst_knots=',num2str(N),'run on_',date,'_at_',num2str(temp(4:5)))
    
        
    blocobj = sampler();
    
    
    %% evaluate and display log posteriors for first iteration
    
    disp('---------------------')
    disp('Un-normalized log posterior for each chain:')
    
    
    for tj = 1:length(tp)
        
        if ~exist('stepvar_chain')
            
            temp_stepvar_chain{tj}           = new_stepvar_chain/tp(tj);
            
        end
        
        [log_r_bot1(tj) PropSoln{tj}] ...
            = logpostprob(dataArray,ddefn,ddehist,N,priorpars,pars_start,tp(tj));
        
        disp(log_r_bot1(tj))
        
    end
    
    if ~exist('stepvar_chain')
        
        stepvar_chain                       = temp_stepvar_chain;
        
    end
    
    disp('---------------------')
    
    for kj=1:kjk
        
        swap_counter{kj} = zeros(1,length(tp));
        
    end
    
    k                       = 1;
    
    %kjk_stop_fix            = 10;                                  % stop tuning after group
    kjk_stop_fix            = floor(kjk/3);                                  % stop tuning after group
    div_by                  = floor(niter/kjk);
    
    for tj = 1:length(tp)
        
        rate{tj}            = zeros(kjk+1,blocobj.blocnum);                    % keep track of the acceptance rate
        accepts{tj}         = zeros(kjk+1,blocobj.blocnum);                    % count when a proposal is accepted
        log_r_bot{tj}       = [log_r_bot1(tj),zeros(1,niter)];
        chainpars{tj}       = [pars_start,zeros(length(pars_start),niter)];    % parameter matrix
        chainSol{tj}        = cell(1,kjk);
        chainSol{tj}{1}     = PropSoln{tj};
        chainTrans{tj}      = cell(1,kjk);
        chainTrans{tj}{1}   = obstransform(PropSoln{tj},chainpars{tj}(5:6,1));
        
    end
    
    
    iter                    = 1;                      % start iteration counter
    a                       = iter;                   % I do this 'a' part because I may stop and restart
    itertime                = 0;                      % keeps track of time per iteration
    
    save(filename)                                    % save everything so far
    
elseif nargin > 3
    
    % load existing file
    tempfilename = filename;
    
    if nargin >= 5
        
        tempcovsname = newcovsfilename;
        
        if nargin == 6
            temprestartiter = restart_at_iter;
        end
        
    end
    
    load(tempfilename);
    
    if nargin == 6
        
        if ~isempty(tempcovsname)
            newcovsfilename = tempcovsname;
            clearvars tempcovsname;
            load(newcovsfilename);
        end
        
        restart_at_iter = temprestartiter;
        clearvars temprestartiter;
        
    elseif nargin == 5
        
        newcovsfilename = tempcovsname;
        clearvars tempcovsname;
        load(newcovsfilename);
        restart_at_iter         = iter+1;
        
    elseif nargin == 4
        
        restart_at_iter         = iter+1;
        
    end
    
    filename                    = tempfilename;
    
    clearvars tempfilename;
    
    
    a                           = restart_at_iter;
    
    if rem(restart_at_iter-1,div_by) == 0
        
        k                        = k+1;
        
    end
    
    for j = 1:length(tp)
        
        rate{j} = [rate{j}(1:k,:);zeros(kjk-k+1,blocobj.blocnum)];                                          % keep track of the acceptance rate
        accepts{j} = [accepts{j}(1:k,:);zeros(kjk-k+1,blocobj.blocnum)];                                         % count when a proposal is accepted
        log_r_bot{j}   = [log_r_bot{j}(1:restart_at_iter+1),zeros(1,niter-restart_at_iter)];                         % log denominator matrix
        chainpars{j} = [chainpars{j}(:,1:restart_at_iter+1),zeros(size(chainpars{j},1),niter-restart_at_iter)];           % parameter matrix
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN THE MCMC CHAIN

for cc = 1:length(tp)
    
    tempchainpars{cc}           = [];
    templog_r_bot{cc}           = [];
    tempaccepts{cc}             = [];
    %tempswap_chains(cc)         = [];
    
end

for iter = a:niter
    
    tic                             % start the iteration timer
    
    rrr                             = iter+1;
    
    %% Propose to swap between chains.
    chain1                          = ceil(rand*(length(tp)+2));
    chain2                          = ceil(rand*(length(tp)+2));
    chains_swap                     = sort([chain1,chain2]);
    
    if chains_swap(1) == chains_swap(2) || chains_swap(1) < 1 || chains_swap(2) > length(chainpars)
        
        chains_swap                 = [0,0];
        
    else
        
        prop_pars{1}                = chainpars{chains_swap(2)}(:,iter);
        prop_pars{2}                = chainpars{chains_swap(1)}(:,iter);
        
        prop_soln{1}                = PropSoln{chains_swap(2)};
        prop_soln{2}                = PropSoln{chains_swap(1)};
        
    end
    
    rr(iter,1:2)                    = chains_swap;
    
    % if there is a swap use these:
    
    if chains_swap(1) > 0
        
        [t_newpars,t_newprobsol,lrb,s_count] = swap_pars_jakstat(...
            {chainpars{chains_swap(1)}(:,iter),chainpars{chains_swap(2)}(:,iter)},...
            {PropSoln{chains_swap(1)},PropSoln{chains_swap(2)}},...
            [log_r_bot{chains_swap(1)}(iter),log_r_bot{chains_swap(2)}(iter)],...
            tp(chains_swap),...
            0,...
            prop_pars,prop_soln,...
            dataArray,priorpars,ddefn,ddehist,N);
        
        
        chainpars{chains_swap(1)}(:,rrr)        = t_newpars{1};
        chainpars{chains_swap(2)}(:,rrr)        = t_newpars{2};
        
        PropSoln{chains_swap(1)}                = t_newprobsol{1};
        PropSoln{chains_swap(2)}                = t_newprobsol{2};
        
        log_r_bot{chains_swap(1)}(rrr)          = lrb(1);
        log_r_bot{chains_swap(2)}(rrr)          = lrb(2);
        
        swap_counter{k}([chains_swap(1),chains_swap(2)])...
            = swap_counter{k}([chains_swap(1),chains_swap(2)]) + s_count;
        
        rr(iter,3)                              = s_count;
        
        
        % if a swap is proposed
        parfor tc = 1:length(tp)                        %non-swap chains only
            
            if tc ~= chains_swap(1) && tc ~= chains_swap(2) && tc < length(tp) + 1
                
                [tempchainpars{tc} tempPropSoln{tc} templog_r_bot{tc} tempaccepts{tc}] ...
                    = sampler(chainpars{tc}(:,iter),PropSoln{tc},log_r_bot{tc}(iter),...
                    stepvar_chain{tc},tp(tc),accepts{tc}(k,:),dataArray,ddefn,ddehist,N,priorpars);
                
            end
            
        end
        
        for tc = 1:length(tp)
            
            if tc ~= chains_swap(1) && tc ~= chains_swap(2) && tc < length(tp) + 1
                
                chainpars{tc}(:,rrr)                = tempchainpars{tc};
                PropSoln{tc}                        = tempPropSoln{tc};
                log_r_bot{tc}(rrr)              	= templog_r_bot{tc};
                accepts{tc}(k,:)                    = tempaccepts{tc};
                
            end
            
        end
        
    else    % If no swap is proposed, then use these:
        
        parfor tc = 1:length(tp) % non-swap chains only
            
            if tc ~= chains_swap(1) && tc ~= chains_swap(2) && tc < length(tp) + 1
                
                [tempchainpars{tc} tempPropSoln{tc} templog_r_bot{tc} tempaccepts{tc}] ...
                    = sampler(chainpars{tc}(:,iter),PropSoln{tc},log_r_bot{tc}(iter),...
                    stepvar_chain{tc},tp(tc),accepts{tc}(k,:),dataArray,ddefn,ddehist,N,priorpars);
                
            end
            
        end
        
        for tc = 1:length(tp)
            %if(tc~=chains_swap(1)&&tc~=chains_swap(2)&& tc<length(tp)+1)
            
            if tc < length(tp)+1
                
                chainpars{tc}(:,rrr)                    = tempchainpars{tc};
                PropSoln{tc}                            = tempPropSoln{tc};
                log_r_bot{tc}(rrr)                      = templog_r_bot{tc};
                accepts{tc}(k,:)                        = tempaccepts{tc};
                
            end
            
        end
        
    end
    
    
    if rem(iter,div_by) == 0
              
        for mm = 1:length(tp)
            
            chainSol{mm}{iter/(div_by)+1}                = PropSoln{mm};
            chainTrans{mm}{iter/(div_by)+1}              = obstransform(PropSoln{mm},chainpars{mm}(5:6,rrr));
            
        end
        
        disp('-----------------')
        disp(iter)
        disp(itertime/iter)
        disp('-----------------')
        
        for mm = 1:length(tp)
            
            rate{mm}(k,:) = ( 1 + accepts{mm}(k,:))./...         % add one accept and one non-accept before computing the ratio.
                (div_by - sum(sum(rr((k-1)*div_by + 1:k*div_by,1:2) == mm)) + 2);
            disp(rate{mm}(k,:))
            
            for aa = 1:blocobj.blocnum
                
                if (k <= kjk_stop_fix) && (rate{mm}(k,aa) < blocobj.bloctarget{aa}(1) || rate{mm}(k,aa) > blocobj.bloctarget{aa}(2))  % Note that I stop making adjustments by kjk_stop_fix iterations
                    
                    stepvar_chain{mm}(blocobj.blocind{aa},blocobj.blocind{aa}) = stepvar_chain{mm}(blocobj.blocind{aa},blocobj.blocind{aa})*rate{mm}(k,aa)/(blocobj.bloctarget{aa}(1) + range(blocobj.bloctarget{aa})/2);
                    
                end
                % target acceptance rates are 44% for one dimension decaying down to 23% for 5 and more dimensions
                
            end
            
        end
        
        adjustcov = floor(kjk_stop_fix/10);
        
        if rem(iter,div_by*adjustcov) == 0 && k <= kjk_stop_fix
             
            try
            
            def_by              = 10;
            span                = (iter-div_by*adjustcov+1):iter;
            
            for ll = 1:length(tp)
                chaincorr{ll}                       = corr(chainpars{ll}(:,span)');
                chaincorr{ll}(isnan(chaincorr{ll})) = 0;   
            end
            
            for ll=1:length(tp)
                
                new_stepcor_chain{ll}     = zeros(npars,npars);
                
                for rr = 1:npars
                    for cc = 1:npars
                        if (rr==cc) 
                            new_stepcor_chain{ll}(rr,cc) = 1;
                        elseif (abs(chaincorr{ll}(rr,cc))<0.6) 
                            new_stepcor_chain{ll}(rr,cc)=0;
                        else
                            new_stepcor_chain{ll}(rr,cc) = sign(chaincorr{ll}(rr,cc))*(abs(chaincorr{ll}(rr,cc))-0.4);
                        end
                    end
                end
            end
            
            
            for ll = 1:length(tp)
                
                pvars{ll} = diag(cov(chainpars{ll}(:,span)'))/def_by;
                
                for rr = 1:npars
                    for cc = 1:npars
                        temp_stepvar_chain{ll}(rr,cc) = new_stepcor_chain{ll}(rr,cc)*sqrt(pvars{ll}(rr)*pvars{ll}(cc));
                    end
                end
                
                tempvar = diag(temp_stepvar_chain{ll});
                
                for dd = 1:npars
                    
                    if tempvar(dd)<=0
                        temp_stepvar_chain{ll}(dd,dd) = stepvar_chain{ll}(dd,dd);
                    end
                    
                end
                
                for bloc = 1:blocobj.blocnum
                    
                    [~,posdef] = chol(temp_stepvar_chain{ll}(blocobj.blocind{bloc},blocobj.blocind{bloc}));
                    
                    if posdef ~= 0 || ~isnan(sum(sum(temp_stepvar_chain{ll}(blocobj.blocind{bloc},blocobj.blocind{bloc}))))
                        stepvar_chain{ll}(blocobj.blocind{bloc},blocobj.blocind{bloc}) = temp_stepvar_chain{ll}(blocobj.blocind{bloc},blocobj.blocind{bloc});
                    end
                    
                end
            end
            
            %disp(stepvar_chain{ll})
            
            catch
                
                %stepvar_chain{ll}
                disp('Error in covariance update')
            
            end
            
        end
        
     save(filename) 
     k                           = k+1;
          
    end
    
    itertime                        = itertime+toc;
    
end

save(filename)

