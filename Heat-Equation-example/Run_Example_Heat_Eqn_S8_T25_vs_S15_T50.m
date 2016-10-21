NumOfIterations = 30; % Increase this to make more accurate

% Run PODES for the heat equation multiple times
for i = 1:NumOfIterations;
    disp(['Running Heat Equation: Iteration ' num2str(i) ' of ' num2str(NumOfIterations)])
    FINAL_PODES_Heat_PDE_Example_Full_Mesh_S8_T25_Coarse
    FINAL_PODES_Heat_PDE_Example_Full_Mesh_S15_T50_Coarse
    %FINAL_PODES_Heat_PDE_Example_Full_Mesh_S29_T100_Coarse
end

% Move to the results folder
cd('./Results')

% Calculate summary statistics of the outputs
Load_Spatial_8_Time_25_Results
Load_Spatial_15_Time_50_Results
%Load_Spatial_29_Time_100_Results

% Plot the summary statistics of the outputs
Plot_Heat_Map_S8_T25_vs_S15_T50
%Plot_Heat_Map_S15_T50_vs_S29_T100

cd('./..')