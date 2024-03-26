% clc; clearvars; close all;

addpath('PyColormap4Matlab')  
% http://www.phy.ohio.edu/~hadizade/blog_files/PyColormap4Matlab.html
% # set the Python path
path_python = '/Users/ywu/opt/anaconda3/envs/tf-env/bin/python';
% # importing Greys colormap from Matplotlib
cl_RdYlBu_r = getPyPlot_cMap('RdYlBu_r', [], [], path_python);
cl_RdBu_r   = getPyPlot_cMap('RdBu_r', [], [], path_python);
cl_RdGy_r   = getPyPlot_cMap('RdGy_r', [], [], path_python);

cmap1 = linspecer;
cmap2 = cl_RdYlBu_r;
cmap3 = cl_RdBu_r;
cmap4 = cl_RdGy_r;
 
%% Define the contour

T = importdata('SputteringCompoMapNormalised.dat');
T_SSS = readtable('MultiTaskModel_v6_NiFeTiMoCr_Retrain_wt_pct_ML_mc_shared_relu.xlsx',"VariableNamingRule","preserve");

T_ML = readtable('MultiTaskModel_v6_NiFeTiMoCr_Retrain_wt_pct_ML_mc_shared_relu.xlsx',"VariableNamingRule","preserve");
Z_H1_ML_mean = T_ML.H1_new_pred_KFold_mean;
Z_H1_ML_std = T_ML.H1_new_pred_KFold_std;

Z_C2_ML_mean = T_ML.C2_new_pred_KFold_mean;
Z_C2_ML_std = T_ML.C2_new_pred_KFold_std;

Z_SSS   = T_SSS.sigma_SSS;

Gmin_eq_FCC_NaN = T_ML.Gmin_eq_FCC;
Gmin_eq_FCC_NaN(Gmin_eq_FCC_NaN==0)=nan;
 
% display(T_ML)

tk_delta = 298.15;
outputname = ['NiFeTiMoCr_', num2str(tk_delta),'K'];
group = [ ...
    {'Ni'}, {'Fe'}, {'Ti'}, {'Mo'}, {'Cr'}; ...
    ];
i = 1;

%% ML_H1_mean
  
figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, Z_H1_ML_mean,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap2)
% caxis([4.8 6.2]);
c = colorbar;
% title('Vickers Hardness','Interpreter', 'latex');
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Vickers hardness by neural network';
c.Label.FontSize = 16;

minZ = floor(min(Z_H1_ML_mean)/100)*100; % Adjust to your data
maxZ = ceil(max(Z_H1_ML_mean)/100)*100;  % Adjust to your data
c.Ticks = minZ:100:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_Retrain_NNH_Hardness_', num2str(tk_delta),'K', '.pdf'],'pdf');


%% ML_H1_mean_FCConly
  
figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, Z_H1_ML_mean.*Gmin_eq_FCC_NaN,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap2)
% caxis([4.8 6.2]);
c = colorbar;
% title('Vickers Hardness','Interpreter', 'latex');
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Vickers hardness by neural network';
c.Label.FontSize = 16;

minZ = floor(min(Z_H1_ML_mean.*Gmin_eq_FCC_NaN)/100)*100; % Adjust to your data
maxZ = ceil(max(Z_H1_ML_mean.*Gmin_eq_FCC_NaN)/100)*100;  % Adjust to your data
c.Ticks = minZ:100:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_Retrain_NNH_Hardness_FCConly_', num2str(tk_delta),'K', '.pdf'],'pdf');


%% ML_H1_std
  
figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, Z_H1_ML_std.*Gmin_eq_FCC_NaN,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap4)
% caxis([4.8 6.2]);
c = colorbar;
% title('Vickers Hardness','Interpreter', 'latex');
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Vickers hardness: Uncertainty';
c.Label.FontSize = 16;

minZ = floor(min(Z_H1_ML_std.*Gmin_eq_FCC_NaN)/100)*100; % Adjust to your data
maxZ = ceil(max(Z_H1_ML_std.*Gmin_eq_FCC_NaN)/100)*100;  % Adjust to your data
c.Ticks = minZ:25:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_Retrain_NNH_Hardness_Uncertainty_', num2str(tk_delta),'K', '.pdf'],'pdf');

 
%% ML_H1_UpperConfidenceBound
  
figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, Z_H1_ML_mean+Z_H1_ML_std,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap4)
% caxis([4.8 6.2]);
c = colorbar;
% title('Vickers Hardness','Interpreter', 'latex');
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Vickers hardness: Upper Confidence Bound';
c.Label.FontSize = 16;

minZ = floor(min(Z_H1_ML_mean+Z_H1_ML_std)/100)*100; % Adjust to your data
maxZ = ceil(max(Z_H1_ML_mean+Z_H1_ML_std)/100)*100;  % Adjust to your data
c.Ticks = minZ:100:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_Retrain_NNH_Hardness_UpperConfidenceBound_', num2str(tk_delta),'K', '.pdf'],'pdf');


%% ML_H1_UpperConfidenceBound_FCConly
  
figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, (Z_H1_ML_mean+Z_H1_ML_std).*Gmin_eq_FCC_NaN,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap4)
% caxis([4.8 6.2]);
c = colorbar;
% title('Vickers Hardness','Interpreter', 'latex');
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Vickers hardness: Upper Confidence Bound';
c.Label.FontSize = 16;

minZ = floor(min(Z_H1_ML_mean+Z_H1_ML_std).*Gmin_eq_FCC_NaN/100)*100; % Adjust to your data
maxZ = ceil(max(Z_H1_ML_mean+Z_H1_ML_std).*Gmin_eq_FCC_NaN/100)*100;  % Adjust to your data
c.Ticks = minZ:100:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_Retrain_NNH_Hardness_UpperConfidenceBound_FCConly_', num2str(tk_delta),'K', '.pdf'],'pdf');



%% ML_C2_mean
  
figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, Z_C2_ML_mean,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap3)
% caxis([4.8 6.2]);
c = colorbar;
% title('Vickers Hardness','Interpreter', 'latex');
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Pitting potential (mV) by neural network';
c.Label.FontSize = 16;

minZ = floor(min(Z_C2_ML_mean)/100)*100; % Adjust to your data
maxZ = ceil(max(Z_C2_ML_mean)/100)*100;  % Adjust to your data
c.Ticks = minZ:100:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_Retrain_NNC_Corrosion_', num2str(tk_delta),'K', '.pdf'],'pdf');




%% ML_C2_mean_FCConly

figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, Z_C2_ML_mean.*Gmin_eq_FCC_NaN,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap3)
% caxis([4.8 6.2]);
c = colorbar;
% title('Vickers Hardness','Interpreter', 'latex');
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')

c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Pitting potential (mV) by neural network';
c.Label.FontSize = 16;

minZ = floor(min(Z_C2_ML_mean.*Gmin_eq_FCC_NaN)/100)*100; % Adjust to your data
maxZ = ceil(max(Z_C2_ML_mean.*Gmin_eq_FCC_NaN)/100)*100;  % Adjust to your data
c.Ticks = minZ:100:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_Retrain_NNC_Corrosion_FCConly_', num2str(tk_delta),'K', '.pdf'],'pdf');


%% ML_C2_std

figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, Z_C2_ML_std.*Gmin_eq_FCC_NaN,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap4)
% caxis([4.8 6.2]);
c = colorbar;
% title('Vickers Hardness','Interpreter', 'latex');
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Pitting potential: Uncertainty';
c.Label.FontSize = 16;

minZ = floor(min(Z_C2_ML_std.*Gmin_eq_FCC_NaN)/100)*100; % Adjust to your data
maxZ = ceil(max(Z_C2_ML_std.*Gmin_eq_FCC_NaN)/100)*100;  % Adjust to your data
c.Ticks = minZ:50:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_Retrain_NNC_Corrosion_Uncertainty_', num2str(tk_delta),'K', '.pdf'],'pdf');


%% ML_C2_UpperConfidenceBound
  
figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, Z_C2_ML_mean+Z_C2_ML_std,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap4)
% caxis([4.8 6.2]);
c = colorbar;
% title('Vickers Hardness','Interpreter', 'latex');
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Pitting potential: Upper Confidence Bound';
c.Label.FontSize = 16;

minZ = floor(min(Z_C2_ML_mean+Z_C2_ML_std)/100)*100; % Adjust to your data
maxZ = ceil(max(Z_C2_ML_mean+Z_C2_ML_std)/100)*100;  % Adjust to your data
c.Ticks = minZ:100:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_Retrain_NNC_Corrosion_UpperConfidenceBound_', num2str(tk_delta),'K', '.pdf'],'pdf');



%% ML_C2_UpperConfidenceBound_FCConly
  
figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, (Z_C2_ML_mean+Z_C2_ML_std).*Gmin_eq_FCC_NaN,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap4)
% caxis([4.8 6.2]);
c = colorbar;
% title('Vickers Hardness','Interpreter', 'latex');
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Pitting potential: Upper Confidence Bound';
c.Label.FontSize = 16;

minZ = floor(min(Z_C2_ML_mean+Z_C2_ML_std).*Gmin_eq_FCC_NaN/100)*100; % Adjust to your data
maxZ = ceil(max(Z_C2_ML_mean+Z_C2_ML_std).*Gmin_eq_FCC_NaN/100)*100;  % Adjust to your data
c.Ticks = minZ:100:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_Retrain_NNC_Corrosion_UpperConfidenceBound_FCConly_', num2str(tk_delta),'K', '.pdf'],'pdf');



%% SSS
  
figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, Z_SSS,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap1)
% caxis([4.8 6.2]);
c = colorbar;
% c.YTick = [200 250 300 350 400 450];
% title('Vickers Hardness','Interpreter', 'latex')
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Solid solution strengthening (MPa)';
c.Label.FontSize = 16;

minZ = floor(min(Z_SSS)/100)*100; % Adjust to your data
maxZ = ceil(max(Z_SSS)/100)*100;  % Adjust to your data
c.Ticks = minZ:50:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_SolidSolutionStrengthening_', num2str(tk_delta),'K', '.pdf'],'pdf');



  %% SSS_FCConly

figure;
set(gcf,'units','points','position',[0,0,600,250]);

hold on
scatter(T(:,1), T(:,2), 300, Z_SSS.*Gmin_eq_FCC_NaN,"filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
scatter(T(:,1).*Gmin_eq_FCC_NaN, T(:,2).*Gmin_eq_FCC_NaN, 300, [0 0 0], 'LineWidth', 2);

colormap(gca,cmap1)
% caxis([4.8 6.2]);
c = colorbar;
% c.YTick = [200 250 300 350 400 450];
% title('Vickers Hardness','Interpreter', 'latex')
xlabel('position x', 'FontSize',14);
ylabel('position y', 'FontSize',14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.7 0.7 0.7], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')
 
c.FontSize = 16;
c.Location = 'eastoutside';
c.Label.String = 'Solid solution strengthening (MPa)';
c.Label.FontSize = 16;

minZ = floor(min(Z_SSS.*Gmin_eq_FCC_NaN)/100)*100; % Adjust to your data
maxZ = ceil(max(Z_SSS.*Gmin_eq_FCC_NaN)/100)*100;  % Adjust to your data
c.Ticks = minZ:100:maxZ;

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
    '-', char(group(i,4)), '-', char(group(i,5)),'_SolidSolutionStrengthening_FCConly_', num2str(tk_delta),'K', '.pdf'],'pdf');




