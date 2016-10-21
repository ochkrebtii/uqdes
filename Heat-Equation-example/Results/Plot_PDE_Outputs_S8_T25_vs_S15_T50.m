function [] = Plot_PDE_Outputs_S8_T25_vs_S15_T50(TimeIndex)

load('Results_PODES_PDE_Spatial_8_Time_25.mat');
load('Results_PODES_PDE_Spatial_15_Time_50.mat');




BoundedLine(0:1/7:1, Mean_Spatial_8_Time_25_Mean(:, TimeIndex+1), Mean_Spatial_8_Time_25_Std(:, TimeIndex+1), '-b');
hold on
BoundedLine(0:1/14:1, Mean_Spatial_15_Time_50_Mean(1:end, 2*TimeIndex+1), Mean_Spatial_15_Time_50_Std(1:end, 2*TimeIndex+1), '-r');



h = errorbar(0:1/7:1, Mean_Spatial_8_Time_25_Mean(:, TimeIndex+1), Mean_Spatial_8_Time_25_Std(:, TimeIndex+1));
set(h,'LineWidth',1.2);
hold on;
h = errorbar(0:1/14:1, Mean_Spatial_15_Time_50_Mean(1:end, 2*TimeIndex+1), Mean_Spatial_15_Time_50_Std(1:end, 2*TimeIndex+1));
set(h,'Color','r');
set(h,'LineWidth',1.2);



% Now plot true solution at the given TimeIndex

% Define the grid points that we want to compare with PODES
kappa = 1;
nt = 25;
dt = 0.25/nt; % Total time = 0.25

nx = 14;
dx = 1/nx; % Total width = 1

% Now refine grid for numerical solver accuracy
Subt = 500;
nt = nt*Subt;
dt = dt/Subt;

Subx = 5;
nx = nx*Subx + 1;
dx = dx/Subx;

% Set initial conditions and spatial grid
T  = sin(0:(pi)/(nx-1):(pi));
x  = 0:dx:dx*(nx-1); %   Grid


time    =   0;
for n=2:nt  % Timestep loop

    % Compute new temperature 
    Tnew   =    zeros(1,nx);
    for i=2:nx-1
        Tnew(i) =  T(i) + (kappa * dt / dx^2) * (T(i+1) - 2*T(i) + T(i-1));
    end

    % Set boundary conditions
    Tnew(1)     =   T(1);
    Tnew(nx)    =   T(nx);
    
    % Update temperature and time
    T           =   Tnew;
    time        =   time+dt;
    
    if mod(n,Subt) == 0
        
        if (n/Subt) == TimeIndex
            h = plot(x(1:Subx:end), Tnew(1:Subx:end), 'g');
            set(h,'LineWidth',1.5);
        end
        
    end
    
end


end