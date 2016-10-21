function [ Results_Mean, Results_Std ] = Load_Spatial_15_Time_50_Results()

SpatialPoints = 15;
TimePoints    = 50;


Files = dir('PODES_PDE_Spatial_15_Time_50_*');

Results = zeros(SpatialPoints*TimePoints,length(Files));

% Load all results
for i = 1:length(Files)
    Result = load(Files(i).name);
    
    Results(:, i) = Result.GP_I_Mean;
end
    

% Get mean and std of results and rescale
Results_Mean = mean(Results, 2);
Results_Std  = std(Results')';


Results_Mean = reshape(Results_Mean, SpatialPoints, TimePoints);
Results_Std  = reshape(Results_Std, SpatialPoints, TimePoints);

Mean_Spatial_15_Time_50_Mean = Results_Mean;
Mean_Spatial_15_Time_50_Std  = Results_Std;

FileName = ['Results_PODES_PDE_Spatial_' num2str(SpatialPoints) '_Time_' num2str(TimePoints)];

save(FileName, 'Mean_Spatial_15_Time_50_Mean', 'Mean_Spatial_15_Time_50_Std');


end

