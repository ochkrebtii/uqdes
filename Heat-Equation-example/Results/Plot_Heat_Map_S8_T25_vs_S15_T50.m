function [ output_args ] = Plot_Heat_Map_S8_T25_vs_S15_T50( )

load('Results_PODES_PDE_Spatial_8_Time_25.mat');
load('Results_PODES_PDE_Spatial_15_Time_50.mat');


subplot(1,3,1)
image([0:1/14:1], [0:0.25/49:0.25], Mean_Spatial_15_Time_50_Mean','CDataMapping', 'scaled')

% Plot full version

subplot(3,3,2)
Plot_PDE_Outputs_S8_T25_vs_S15_T50(2)
axis([-0.01 1.01 0 0.85])

subplot(3,3,2+3)
Plot_PDE_Outputs_S8_T25_vs_S15_T50(12)
axis([-0.01 1.01 0 0.35])
axis([-0.01 1.01 0 0.85])

subplot(3,3,2+6)
Plot_PDE_Outputs_S8_T25_vs_S15_T50(22)
axis([-0.01 1.01 0 0.18])
axis([-0.01 1.01 0 0.85])

% Plot zoomed in version

subplot(3,3,3)
Plot_PDE_Outputs_S8_T25_vs_S15_T50(2)
axis([0.4 0.6 0.74 0.88])

subplot(3,3,3+3)
Plot_PDE_Outputs_S8_T25_vs_S15_T50(12)
axis([0.4 0.6 0.22 0.36])

subplot(3,3,3+6)
Plot_PDE_Outputs_S8_T25_vs_S15_T50(22)
axis([0.4 0.6 0.03 0.17])


end

