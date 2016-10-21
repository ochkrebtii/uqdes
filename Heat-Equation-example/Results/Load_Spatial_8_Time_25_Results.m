function [ Results_Mean, Results_Std ] = Load_Spatial_8_Time_25_Results()

SpatialPoints = 8;
TimePoints    = 25;


Files = dir('PODES_PDE_Spatial_8_Time_25_*');

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


Mean_Spatial_8_Time_25_Mean = Results_Mean;
Mean_Spatial_8_Time_25_Std  = Results_Std;

FileName = ['Results_PODES_PDE_Spatial_' num2str(SpatialPoints) '_Time_' num2str(TimePoints)];

save(FileName, 'Mean_Spatial_8_Time_25_Mean', 'Mean_Spatial_8_Time_25_Std');

end

