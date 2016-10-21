function [] = FINAL_PODES_Heat_PDE_Example_Full_Mesh_S15_T50_Coarse()

addpath(genpath(pwd)) % Add all subfolders in this folder

% Spatial discretisation
x_start               = 0;
x_end                 = 1;
NumOfSpatialPoints    = 15;
int_step_x            = (x_end-x_start)/(NumOfSpatialPoints-1);
Spatial               = x_start:int_step_x:x_end; % Spatial integration points


NumOfTimePoints       = 50;
int_step              = 0.25/NumOfTimePoints;
Time                  = 0:int_step:((NumOfTimePoints-1)*int_step);


Initial_Conditions    = sin(0:(pi)/(NumOfSpatialPoints-1):(pi)); % Column vector


StateData             = zeros(NumOfSpatialPoints, NumOfTimePoints);
State_Noise_Var       = zeros(NumOfSpatialPoints, NumOfTimePoints);
Deriv_Data            = zeros(NumOfSpatialPoints, NumOfTimePoints);
Deriv_Noise_Var       = zeros(NumOfSpatialPoints, NumOfTimePoints);


% Observations of the function at initial conditions
StateData(:, 1)       = Initial_Conditions;
State_Noise_Var(:,1)  = ones(NumOfSpatialPoints, 1)*1e-8;
% Heat equation
Deriv_Data(:, 1)      = -sin(Spatial*pi)*pi^2; % Analytic derivative of time is second derivative wrt space
Deriv_Noise_Var(:, 1) = ones(NumOfSpatialPoints, 1)*1e-8;




%% PDE PODES %%

% Time Covariance Functions
HP_Time  = [1.5*int_step, 0.0001]; % uniform bounded cov
CC_Time  = RR1d(Time, Time, HP_Time(1), Time(1), Time(end))./HP_Time(2);
CI_Time  = RQ1d(Time, Time, HP_Time(1), Time(1), Time(end))./HP_Time(2);
II_Time  = QQ1d(Time, Time, HP_Time(1), Time(1), Time(end))./HP_Time(2);

% Spatial Covariance Functions
HP_Space   = [2, 6*int_step_x]; % squared exponential cov
CC_Space   = RBF_CC(HP_Space, Spatial, Spatial);
CD2_Space  = RBF_CD2(HP_Space, Spatial, Spatial);
D2D2_Space = RBF_D2D2(HP_Space, Spatial, Spatial);

% Set up initial noise covariance
NoiseCov = ones(NumOfSpatialPoints, 1)*1e-8;

% Set spatial indices
Idx = 1:NumOfSpatialPoints;


for i = 2:NumOfTimePoints
    
    %disp(['Time Step: ' num2str(i)])
    
    %%% Calculate derivative estimates at current time point i %%%
    if i == 2
        % Calculate exactly
        Deriv_Data(:,2) = -sin(Spatial*pi)*pi^2; % Exact
    else
        
        % Estimate using GP predictive distribution
        %
        A_1_A = kron(CC_Time(1:i-1, 1:i-1), CC_Space(Idx, Idx));
        A_1_C = kron( CI_Time(1:i-1,2:i), CC_Space(Idx, [1 NumOfSpatialPoints]) );
        A_1_B = kron( II_Time(2:i,2:i), CC_Space([1 NumOfSpatialPoints], [1 NumOfSpatialPoints]) );
        A_1 = [A_1_A A_1_C; A_1_C' A_1_B];

        A_2_A = kron(CI_Time(1:i-1, i), CD2_Space(Idx, Idx));
        A_2_B = kron( II_Time(2:i, i), CD2_Space([1 NumOfSpatialPoints], Idx) );
        A_2 = [A_2_A; A_2_B];

        A_3 = D2D2_Space(Idx, Idx).*II_Time(i,i);
        

        A = A_1;
        C = A_2;
        B = A_3;

        Obs_Values = Deriv_Data(Idx,1:i-1); % Get all
        Obs_Values = Obs_Values(:)'; % Concatenate
        Int_Obs = zeros(1, 2*length(2:i));
        Int_NoiseCov = ones(1, 2*length(2:i))*1e-6;

        % Modelling state as IC + GP, therefore derivative data is d^2/dx^2(IC + GP)
        GP_D2_Mean = (C'*( ( A + diag(sparse([NoiseCov' Int_NoiseCov])) )\[Obs_Values Int_Obs]' ));
        GP_D2_Var  = B - C'/(A+ diag(sparse([NoiseCov' Int_NoiseCov])) )*C;
        
        %Deriv_Data(:,i) = Deriv_Data(:, 1) + GP_D2_Mean;
        Deriv_Data(:,i) = Deriv_Data(:, 1) + GP_D2_Mean + ( randn(1, NumOfSpatialPoints)*chol(GP_D2_Var+eye(NumOfSpatialPoints)*1e-6) )';
        Deriv_Data([1 NumOfSpatialPoints], i) = [0 0];
        
        %}
    end
    
    % Predict variance of derivative observations at next time point - mean and variance
    
    A_1_A = kron(CC_Time(1:i-1, 1:i-1), CC_Space(Idx, Idx));
    A_1_C = kron( CI_Time(1:i-1,2:i), CC_Space(Idx, [1 NumOfSpatialPoints]) );
    A_1_B = kron( II_Time(2:i,2:i), CC_Space([1 NumOfSpatialPoints], [1 NumOfSpatialPoints]) );
    A_1 = [A_1_A A_1_C; A_1_C' A_1_B];
    
    A_2_A = kron(CC_Time(1:i-1, i), CC_Space(Idx, Idx));
    A_2_B = kron( CI_Time(i,2:i), CC_Space(Idx, [1 NumOfSpatialPoints]) )';
    A_2 = [A_2_A; A_2_B];
    
    A_3 = CC_Space(Idx, Idx).*CC_Time(i,i);
    
    
    A = A_1;
    C = A_2;
    B = A_3;
    
    Obs_Values = Deriv_Data(Idx,1:i-1); % Get all
    Obs_Values = Obs_Values(:)'; % Concatenate
    Int_Obs = zeros(1, 2*length(2:i));
    Int_NoiseCov = ones(1, 2*length(2:i))*1e-6;
    
    %GP_C_Mean = C'*( ( A + diag(sparse([NoiseCov' Int_NoiseCov])) )\[Obs_Values Int_Obs]' );
    GP_C_Var  = B - C'/(A+ diag(sparse([NoiseCov' Int_NoiseCov])) )*C;
    
    NewNoiseCov = diag(GP_C_Var);
    NoiseCov = [NoiseCov; 0; NewNoiseCov(2:NumOfSpatialPoints-1); 0];
    
    
    
end





A_1_A = kron(CC_Time(1:i, 1:i), CC_Space(Idx, Idx));
A_1_C = kron( CI_Time(1:i,2:NumOfTimePoints), CC_Space(Idx, [1 NumOfSpatialPoints]) );
A_1_B = kron( II_Time(2:NumOfTimePoints,2:NumOfTimePoints), CC_Space([1 NumOfSpatialPoints], [1 NumOfSpatialPoints]) );
A_1 = [A_1_A A_1_C; A_1_C' A_1_B];

A_2_A = kron(CI_Time(1:i, 1:NumOfTimePoints), CC_Space(Idx, Idx));
A_2_B = kron( II_Time(1:NumOfTimePoints,2:NumOfTimePoints), CC_Space(Idx, [1 NumOfSpatialPoints]) )';
A_2 = [A_2_A; A_2_B];

A_3 = kron( II_Time(1:NumOfTimePoints,1:NumOfTimePoints), CC_Space(Idx, Idx) );


A = A_1;
C = A_2;
B = A_3;

Obs_Values = Deriv_Data(Idx,1:i); % Get all
Obs_Values = Obs_Values(:)'; % Concatenate
Int_Obs = zeros(1, 2*length(2:NumOfTimePoints));
Int_NoiseCov = ones(1, 2*length(2:NumOfTimePoints))*1e-6;

GP_I_Mean = C'*( ( A + diag(sparse([NoiseCov' Int_NoiseCov])) )\[Obs_Values Int_Obs]' );
GP_I_Var  = B - C'/(A+ diag(sparse([NoiseCov' Int_NoiseCov])) )*C;

GP_I_Mean = reshape(GP_I_Mean, NumOfSpatialPoints, NumOfTimePoints) + repmat(StateData(Idx, 1), 1, NumOfTimePoints);
GP_I_Mean = reshape(GP_I_Mean, NumOfSpatialPoints*NumOfTimePoints, 1);

TimeStamp = ceil(now*100000);
CurrentFileID = ['PODES_PDE_Spatial_' num2str(NumOfSpatialPoints) '_Time_' num2str(NumOfTimePoints) '_' num2str(TimeStamp)];
save(['./Results/' CurrentFileID], 'GP_I_Mean', 'GP_I_Var', 'HP_Time', 'HP_Space');




end