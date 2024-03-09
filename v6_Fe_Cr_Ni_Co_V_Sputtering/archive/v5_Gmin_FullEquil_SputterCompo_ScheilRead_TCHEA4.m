% 20200407
% Calculation of Fe-C-V-Al alloy V content in cementite at 575C

clc,clearvars,close all

% addpath('H:\Matlab Toolbox HEA\FunctionLib');
% % % add_ternary_paths

%% Define the contour
T = importdata('SputteringCompoMapNormalised.dat');
tk_delta = 298.15;
% tk      = (900+273.15) : 50 : (1300+273.15);

% outputname = ['KW130_FeCrNi-MoTi_', num2str(tk_delta),'K'];

group = [ ...
    {'Ni'}, {'Cr'}, {'Co'}, {'V'}, {'Fe'}; ...
    ];

%% Define contour matrix and start calculation

% for i = 1: size(group,1)
i =1;
Z_ScheilFCC = NaN(size(T,1),1);



%% Select elements and retrieve data
sel=[char(group(i,1)), ' ', ...
    char(group(i,2)), ' ', ...
    char(group(i,3)), ' ', ...
    char(group(i,4)), ' ', ...
    char(group(i,5)),  ...
    ];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print out Scheil INPUT FILE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:length(T)

    x1 = T(j,3) * 100;
    x2 = T(j,4) * 100;
    x3 = T(j,5) * 100;
    x4 = T(j,6) * 100;
    x5 = T(j,7) * 100;

    scheil_fname  = ['Scheil_TCHEA4_P',num2str(j),'_', num2str(x1),'Ni-', num2str(x2),'Cr-',...
        num2str(x3),'Co-', num2str(x4),'V-',num2str(x5),'Fe'];
    old = ".";
    new = "_";
    scheil_fname = replace(scheil_fname,old,new);
    scheil_OUTPUT_fname = ['OUTPUT_', scheil_fname, '_FCC_NS_T.exp'];


    if isfile(['NiCrCoVFe_Scheil_TCHEA4/', scheil_OUTPUT_fname])

        ScheilCurve_temp = NumExtn({['NiCrCoVFe_Scheil_TCHEA4/', scheil_OUTPUT_fname]});
        ScheilCurve_temp = unique(ScheilCurve_temp,'rows');

        % ScheilFreeze_temp = ScheilCurve_temp(1,2) - ScheilCurve_temp(end,2);
        % Z_ScheilFreeze(imo,ini) = ScheilFreeze_temp;
        % T_ScheilFreeze = [T_ScheilFreeze; [wni(ini) wmo(imo) ScheilFreeze_temp] ];

        ScheilFCC_temp = ScheilCurve_temp(end,1);
        Z_ScheilFCC(j) = ScheilFCC_temp;

    end
end

%% Plotting - all

cmap = linspecer;

figure;
set(gcf,'units','points','position',[0,0,800,300]);

subplot(1, 2, 1)
hold on
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5])
scatter(T(:,1), T(:,2), 300, Z_ScheilFCC, "filled")
colormap(gca,cmap)
caxis([0.5 1]);
colorbar;
title('TCHEA4: Scheil solidification FCC fraction','Interpreter', 'latex');
xlabel('position x');
ylabel('position y');
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',16); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',16);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, 'horizontal','center', 'vertical','middle', 'color', [0.5 0.5 0.5])

subplot(1, 2, 2)
hold on
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5])
scatter(T(:,1), T(:,2), 300, Z_ScheilFCC, "filled")
colormap(gca,cmap)
caxis([0.90 1]);
colorbar;
title('TCHEA4: Scheil solidification FCC fraction ($>$90$\%$)','Interpreter', 'latex');
xlabel('position x');
ylabel('position y');
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize',16); xtickangle(45);
set(gca,'yscale','lin', 'FontSize',16);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, 'horizontal','center', 'vertical','middle', 'color', [0.5 0.5 0.5])

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
 '-', char(group(i,4)), '-', char(group(i,5)),'_Scheil_TCHEA4_FCC_', num2str(tk_delta),'K', '.png'],'png');




%% Plotting - single

cmap = linspecer;

figure;
set(gcf,'units','points','position',[0,0,400,250]);

subplot(1, 1, 1)
hold on
scatter(T(:,1), T(:,2), 300, Z_ScheilFCC, "filled")
scatter(T(:,1), T(:,2), 300, [0.5 0.5 0.5])
colormap(gca,cmap)
caxis([0 1]);

c = colorbar;
c.FontSize = 16;
c.Label.String = 'FCC fraction (Scheil simulation)';
c.Label.FontSize = 16;
c.YTick = [0 0.2 0.4 0.6 0.8 1];
c.YTickLabel = {'0', '20%', '40%', '60%', '80%', '100%'};
% c.YTick = [0.9 0.92 0.94 0.96 0.98 1];
% c.YTickLabel = {'<90%', '92%', '94%', '96%', '98%', '100%'};

% title('TCHEA4: Scheil solidification FCC fraction ($>$90$\%$)','Interpreter', 'latex');
xlabel('position x', 'FontSize', 14);
ylabel('position y', 'FontSize', 14);
box on; axis square; xticks(0:10:100); yticks(0:10:100);
set(gca,'xscale','lin', 'FontSize', 14); xtickangle(45);
set(gca,'yscale','lin', 'FontSize', 14);

labels = num2str((1:size(T,1))','%d');
text(T(:,1), T(:,2), labels, ...
    'horizontal','center', 'vertical','middle', ...
    'color', [0.8 0.8 0.8], ...
    'FontSize', 12, ...
    'FontWeight', 'bold')

saveas(gcf, [char(group(i,1)), '-',char(group(i,2)), '-', char(group(i,3)), ...
 '-', char(group(i,4)), '-', char(group(i,5)),'_Scheil_TCHEA4_FCC_', num2str(tk_delta),'K-single', '.pdf'],'pdf');




 