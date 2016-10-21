function [ output_args ] = Plot_Heat_Map_S15_T50_vs_S29_T100( )

load('Results_PODES_PDE_Spatial_15_Time_50.mat')
load('Results_PODES_PDE_Spatial_29_Time_100.mat')

subplot(1,3,1)
image([0:1/14:1], [0:0.25/49:0.25], Mean_Spatial_15_Time_50_Mean','CDataMapping', 'scaled')

% Plot full version

subplot(3,3,2)
Plot_PDE_Outputs_S15_T50_vs_S29_T100(4)
axis([-0.01 1.01 0 0.85])

subplot(3,3,2+3)
Plot_PDE_Outputs_S15_T50_vs_S29_T100(24)
axis([-0.01 1.01 0 0.35])
axis([-0.01 1.01 0 0.85])

subplot(3,3,2+6)
Plot_PDE_Outputs_S15_T50_vs_S29_T100(44)
axis([-0.01 1.01 0 0.18])
axis([-0.01 1.01 0 0.85])

% Plot zoomed in version

subplot(3,3,3)
Plot_PDE_Outputs_S15_T50_vs_S29_T100(4)
axis([0.4 0.6 0.78 0.84])

subplot(3,3,3+3)
Plot_PDE_Outputs_S15_T50_vs_S29_T100(24)
axis([0.4 0.6 0.26 0.32])

subplot(3,3,3+6)
Plot_PDE_Outputs_S15_T50_vs_S29_T100(44)
axis([0.4 0.6 0.07 0.13])


end

