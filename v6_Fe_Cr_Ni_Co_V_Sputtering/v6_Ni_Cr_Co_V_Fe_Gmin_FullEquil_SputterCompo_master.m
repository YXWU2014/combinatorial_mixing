clc,clearvars,close all
addpath('H:\Matlab Toolbox HEA\FunctionLib');

% http://www.phy.ohio.edu/~hadizade/blog_files/PyColormap4Matlab.html
% # set the Python path
addpath('.\PyColormap4Matlab')
% path_python = '/Users/ywu/opt/anaconda3/envs/tf-env/bin/python';
path_python = 'C:\Users\yuxiang.wu\AppData\Local\anaconda3\python.exe';
% # importing Greys colormap from Matplotlib
cl_Blues = getPyPlot_cMap('Blues', [], [], path_python);
cl_RdBu_r  = getPyPlot_cMap('RdBu_r', [], [], path_python);
cl_RdGy_r  = getPyPlot_cMap('RdGy_r', [], [], path_python);

cmap1 = linspecer;
cmap2 = cl_Blues;
cmap3 = cl_RdBu_r;
cmap4 = cl_RdGy_r;

%% Define the input
% PVD_compo = readmatrix('SputteringCompoMapNormalised.dat');
PVD_compo = importdata('SputteringCompoMapNormalised.dat');

[~, filename, ~] = fileparts(mfilename('fullpath'));

% Split filename into parts using '_' as delimiter
filename_parts = strsplit(filename, '_');
ele_A = filename_parts{2};
ele_B = filename_parts{3};
ele_C = filename_parts{4};
ele_D = filename_parts{5};
ele_E = filename_parts{6};

% group = perms([{ele_A} {ele_B} {ele_C} {ele_D} {ele_E}]);
group = [{ele_A} {ele_B} {ele_C} {ele_D} {ele_E}];

compo_table =  table([],[],[],[],[],[],[],[], ...
    'VariableNames', {ele_A, ele_B, ele_C, ele_D, ele_E, 'Gmin_FCC', 'eq_FCC', 'Gmin_eq_FCC'});

count_fcc_fractions = nan(size(group,1),3);
count_ABCDE         = char(nan(size(group,1), length([char(group(1,1)), '-',char(group(1,2)), ...
    '-', char(group(1,3)), '-', char(group(1,4)), '-', char(group(1,5))])));

%%
for permu_i = 1: size(group,1)

    %     % Skip the calculated
    %     if isfile([char(group(permu_i,1)), '-',char(group(permu_i,2)),'-', char(group(permu_i,3)),...
    %             '-', char(group(permu_i,4)), '-', char(group(permu_i,5)), '_Sputter compo contour_', num2str(tk_Gmin),'K_','.jpeg'])
    %         continue;
    %     end

    % permu_i = 1

    %% plotting the sputtering composition contour
    % close all;
    figure(1)
    % caxis_limits = [20 60; 5 45; 0 12; 0 25; 10 70]; % Define caxis limits for each subplot

    for subplot_j = 1:5
        subplot(2, 3, subplot_j)
        hold on
        scatter(PVD_compo(:,1), PVD_compo(:,2), 300, PVD_compo(:, subplot_j+2)*100, "filled")
        scatter(PVD_compo(:,1), PVD_compo(:,2), 300, [0.5 0.5 0.5],'LineWidth',1)
        colormap(gca,'parula')
        c=colorbar;
        title(char(group(permu_i, subplot_j)), 'FontWeight', 'normal');
        xlabel('position x', 'FontSize',14);
        ylabel('position y', 'FontSize',14);

        % caxis(caxis_limits(subplot_j, :)) % Set caxis limits

        box on; axis square;
        xticks(0:10:100); yticks(0:10:100);
        set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
        set(gca,'yscale','lin', 'FontSize',14);
        labels = num2str((1:size(PVD_compo,1))','%d');
        text(PVD_compo(:,1), PVD_compo(:,2), labels, ...
            'horizontal','center', 'vertical','middle', ...
            'color', [0.5 0.5 0.5], ...
            'FontSize', 12, ...
            'FontWeight', 'bold')

        c.FontSize =16;
        c.Label.String = 'at.%';
        c.Label.FontSize = 16;
        c.Label.Rotation = 0; % Makes the label horizontal
        c.Label.Position = [1 max(get(c,'Limits'))*1.12 0]; % Adjusts the label to top
    end

    set(gcf,'units','points','position',[0,0,1200,550]);
    set(gcf,'PaperSize',[60 40]);

    % saveas(gcf, [char(group(permu_i,1)), '-',char(group(permu_i,2)), '-', char(group(permu_i,3)), '-', ...
    %     char(group(permu_i,4)), '-', char(group(permu_i,5)),'_Sputter_CompoContour','.pdf'],'pdf');
    print(gcf,[char(group(permu_i,1)), '-',char(group(permu_i,2)), '-', char(group(permu_i,3)), '-', ...
        char(group(permu_i,4)), '-', char(group(permu_i,5)),'_Sputter_CompoContour','.pdf'], '-dpdf')


    %% Initiate TC system
    tc_init_root;
    tc_open_database('tchea4');
    tc_check_error;

    % Select elements and retrieve data
    sel=[char(group(permu_i,1)), ' ', char(group(permu_i,2)), ' ', ...
        char(group(permu_i,3)), ' ', char(group(permu_i,4)), ' ', char(group(permu_i,5))];
    tc_element_select(sel);
    % phases = 'diamond graphite';
    % tc_phase_reject(phases)
    tc_get_data;

    % Set conditions
    tc_set_condition('n',1);
    tc_set_condition('p',101325);

    %% --- calculate minimum G phase diagram

    % Define the list of phase names
      phase_names = {'BCC_B2', 'DIS_FCC_A1', 'FCC_L12', 'HEUSLER_L21', 'SIGMA', ...
        'LIQUID', 'NI3TA_D0A', 'MU_PHASE', 'C14_LAVES', 'C15_LAVES', 'CHI_A12'};

    tk_Gmin = 473.15;

    % Initialize the output arrays
    n_compositions = length(PVD_compo);
    Z_gm = cell(1, length(phase_names));
    Z_gm_Gmin = cell(1, length(phase_names));

    for k = 1:length(phase_names)
        Z_gm{k} = [PVD_compo(:,[1 2]), NaN(n_compositions,1)];
        Z_gm_Gmin{k} = [PVD_compo(:,[1 2]), NaN(n_compositions,1)];
    end

    % Loop over each composition
    for j = 1:n_compositions
        % Set the composition and temperature conditions
        tc_set_condition(['x(', char(group(permu_i,1)), ')'], PVD_compo(j,3));
        tc_set_condition(['x(', char(group(permu_i,2)), ')'], PVD_compo(j,4));
        tc_set_condition(['x(', char(group(permu_i,3)), ')'], PVD_compo(j,5));
        tc_set_condition(['x(', char(group(permu_i,4)), ')'], PVD_compo(j,6));
        tc_set_condition('T', tk_Gmin);

        % Loop over each phase
        gm_temp = nan(1, length(phase_names));
        for k = 1:length(phase_names)
            % Set the phase status and minimize the Gibbs energy
            tc_set_phase_status('*', 'suspended', 0);
            tc_set_phase_status(phase_names{k}, 'entered', 1);
            tc_set_minimization('off');
            tc_compute_equilibrium();

            % Get the Gibbs energy and update the output array
            [err_code, err_msg] = tc_error();
            if err_code == 0
                gm_temp(k) = tc_get_value(['gm(', phase_names{k}, ')']);
                Z_gm{k}(j, 3) = gm_temp(k);
            else
                fprintf('Error %d: %s\n', err_code, err_msg);
            end
            tc_reset_error();
        end

        % Find the minimum Gibbs energy and update the output array
        [~, idx] = min(gm_temp);
        Z_gm_Gmin{idx}(j, 3) = gm_temp(idx);
    end

    %% plot: calculate minimum G phase diagram

    % Create the figure and loop over the GMs
    figure('units', 'points', 'position', [0, 0, 2000, 800]);
    marker_size = 75;
    
    for phase_names_k = 1:length(phase_names)
        
        % Define the subplot location: gm
        subplot(4, 6, phase_names_k);
        % Scatter plot the data and apply settings
        scatter(PVD_compo(:, 1), PVD_compo(:, 2), marker_size, Z_gm{phase_names_k}(:, 3), 'filled');
        c = colorbar;
        c.Label.String = 'Molar Gibbs energy-gm (J/mol)';
        c.Ruler.Exponent=4;                   % set the desired exponent
        c.Ruler.TickLabelFormat='%0.1f';       % fix up ugly default %g formatting
        % title(['gm(', phase_names{phase_names_k}, ')'], 'Interpreter', 'none');
        titleHandle = title(['gm(', phase_names{phase_names_k}, ')'], 'Interpreter', 'none');
        titlePos = get(titleHandle, 'Position');
        set(titleHandle, 'Position', [titlePos(1), titlePos(2) * 1.05, titlePos(3)]);
        xlabel('position x');
        ylabel('position y');
        box on; axis square;
        xticks(0:10:100); xtickangle(45);yticks(0:10:100);
       
        
        % Define the second subplot location: gm_Gmin
        subplot(4, 6, phase_names_k + 12);
        % Scatter plot the second data and apply settings
        scatter(PVD_compo(:, 1), PVD_compo(:, 2), marker_size, Z_gm_Gmin{phase_names_k}(:, 3), 'filled');
        c = colorbar;
        c.Label.String = 'Molar Gibbs energy-gm (J/mol)';
        caxis([-2.3e4 -1.5e4]) 
        % title(['gm(', phase_names{phase_names_k}, ')_Gmin'], 'Interpreter', 'none');
        titleHandle = title(['gm(', phase_names{phase_names_k}, ')_Gmin'], 'Interpreter', 'none');
        titlePos = get(titleHandle, 'Position');
        set(titleHandle, 'Position', [titlePos(1), titlePos(2) * 1.05, titlePos(3)]);
        
        xlabel('position x'); ylabel('position y');
        box on; axis square;
        xticks(0:10:100); xtickangle(45);yticks(0:10:100);
    end
    
    set(gcf,'PaperSize',[30 20]);

    % Save the figure for the current GM
    print(gcf, [char(group(permu_i,1)), '-', char(group(permu_i,2)), '-', char(group(permu_i,3)), ...
        '-', char(group(permu_i,4)), '-', char(group(permu_i,5)), '_G comparison_', num2str(tk_Gmin), 'K.pdf'], '-dpdf');

    %% full equilibrium phase diagram

    % Z_eq to store if each phase has single-phase range
    Z_eq = cell(1, length(phase_names));
    for phase_names_k = 1:length(phase_names)
        Z_eq{phase_names_k} = [PVD_compo(:,[1 2]), NaN(n_compositions,1)];
    end

    for j = 1: length(PVD_compo)
        tc_set_condition(['x(', char(group(permu_i,1)), ')'], PVD_compo(j,3));
        tc_set_condition(['x(', char(group(permu_i,2)), ')'], PVD_compo(j,4));
        tc_set_condition(['x(', char(group(permu_i,3)), ')'], PVD_compo(j,5));
        tc_set_condition(['x(', char(group(permu_i,4)), ')'], PVD_compo(j,6));
        % tc_set_condition(['x(', char(group(permu_i,5)), ')'], PVD_compo(j,7));

        % looping the temperature
        tk = (900+273.15) : 50 : (1300+273.15);
        np_eq_tk = cell(1, length(phase_names));

        % np_eq_tk to loop temperature range for each phase at one PVD_compo
        for phase_names_k = 1:length(phase_names)
            np_eq_tk{phase_names_k} = NaN(length(tk), 1);
        end

        for i_tk = 1:length(tk) % looping the temperature

            tc_set_condition('T', tk(i_tk));
            tc_set_phase_status('*', 'entered', 1);
            tc_set_minimization('off')

            tc_compute_equilibrium;
            [a,b]=tc_error;
            if a~=0
                s=sprintf(' ERROR %d',a);
                disp(s);
                disp(b);
            else

                % store the np at certain temperature for each phase
                for phase_names_k = 1:length(phase_names)
                    np_eq_tk{phase_names_k}(i_tk)  = tc_get_value(['np(', phase_names{phase_names_k},')']);

                    if strcmp(phase_names{phase_names_k}, 'FCC_L12') % this is phase name for FCC_L12
                        if np_eq_tk{phase_names_k}(i_tk) == 0
                            np_eq_tk{phase_names_k}(i_tk) = tc_get_value(['np(', phase_names{phase_names_k},'#1)']);
                        end
                        if np_eq_tk{phase_names_k}(i_tk) == 0
                            np_eq_tk{phase_names_k}(i_tk) = tc_get_value(['np(', phase_names{phase_names_k},'#2)']);
                        end
                    end
                end

            end
            tc_reset_error;
        end

        % decide if single phase for each phase in the temperature range
        for phase_names_k = 1:length(phase_names)

            np_eq_tk_temp = np_eq_tk{phase_names_k}; % temp variable

            if  ~isnan(np_eq_tk_temp) % if np_eq_tk with number
                np_eq_tk_temp = np_eq_tk_temp(np_eq_tk_temp >= 0.99); % with 99% np to be existing
                Z_eq{phase_names_k}(j,3) = ~isempty(np_eq_tk_temp); % if not empty=1, empty=0
            end
        end
    end

    %% extract the FCC phase data
    Z_gm_Gmin_FCC      = Z_gm_Gmin{find(strcmp(phase_names, 'FCC_L12'))};
    Z_gm_Gmin_FCC(:,3) = ~isnan(Z_gm_Gmin_FCC(:,3));

    Z_eq_FCC           = Z_eq{find(strcmp(phase_names, 'FCC_L12'))};

    Z_gm_Gmin_eq_FCC = [PVD_compo(:,[1 2]), zeros(n_compositions,1)];
    Z_gm_Gmin_eq_FCC(:,3) = Z_gm_Gmin_FCC(:,3) == 1 & Z_eq_FCC(:,3) == 1;

    compo_data_j = table(PVD_compo(:, 3), PVD_compo(:, 4), PVD_compo(:, 5), PVD_compo(:, 6), PVD_compo(:, 7),...
        Z_gm_Gmin_FCC(:, 3), Z_eq_FCC(:, 3), Z_gm_Gmin_eq_FCC(:,3), ...
        'VariableNames', {char(group(permu_i,1)), char(group(permu_i,2)), char(group(permu_i,3)), char(group(permu_i,4)), char(group(permu_i,5)), ...
        'Gmin_FCC', 'eq_FCC', 'Gmin_eq_FCC'});

    compo_data_j = movevars(compo_data_j, ele_B, 'After', ele_A);
    compo_data_j = movevars(compo_data_j, ele_C, 'After', ele_B);
    compo_data_j = movevars(compo_data_j, ele_D, 'After', ele_C);
    compo_data_j = movevars(compo_data_j, ele_E, 'After', ele_D);

    compo_table = vertcat(compo_table, compo_data_j);

    %% plotting: minimum G diagram
    figure(3);
    set(gcf, 'Units', 'points', 'Position', [0, 0, 700, 250]);
    cmap = linspecer(length(phase_names));
    % t = tiledlayout('flow', 'TileSpacing', 'compact');

    % Create scatter plot
    % nexttile
    hold on;
    phase_names_sel = int16.empty;
    for phase_names_k = 1:length(phase_names)
        LineNoNaN = find(~isnan(Z_gm_Gmin{phase_names_k}(:,3)));
        
        if ~isnan(LineNoNaN)
           phase_names_sel(end+1) = phase_names_k;
        end 
        
        h(phase_names_k) = scatter(PVD_compo(LineNoNaN,1), PVD_compo(LineNoNaN,2), 300, cmap(phase_names_k,:), 'filled');
    end
    LineNoNaN_overlap = find(Z_gm_Gmin_eq_FCC(:,3)==1);
    h_ = scatter(PVD_compo(LineNoNaN_overlap,1), PVD_compo(LineNoNaN_overlap,2), 300, [0 0.4470 0.7410], 'LineWidth', 2);

    xlabel('position x', 'FontSize', 14);
    ylabel('position y', 'FontSize', 14);
    box on;
    axis square;
    xticks(0:10:100);
    xtickangle(45);
    yticks(0:10:100);
    set(gca, 'XScale', 'lin', 'YScale', 'lin', 'FontSize', 14);
    hold off;

    % Add labels and legend
    labels = num2str((1:size(PVD_compo,1))','%d');
    text(PVD_compo(:,1), PVD_compo(:,2), labels, ...
        'horizontal','center', 'vertical','middle', ...
        'color', [0.3 0.3 0.3], ...
        'FontSize', 12, ...
    'FontWeight', 'bold')
    legend(h(phase_names_sel), phase_names(phase_names_sel), ...
        'Interpreter', 'none', ...
        'Location', 'bestoutside', ...
       'FontSize', 16);

    % Save plot
    saveas(gcf, [char(group(permu_i,1)), '-', char(group(permu_i,2)), '-', char(group(permu_i,3)), ...
        '-', char(group(permu_i,4)), '-', char(group(permu_i,5)), '_Gmin-FullEquil_Gmin_', num2str(tk_Gmin), 'K', '.pdf'], 'pdf');

    %% --- plotting: full equilibrium phase diagram
    figure(4);
    set(gcf,'units','points','position',[0,0,600,250]);

    % t = tiledlayout('flow','TileSpacing','compact');
    % nexttile

    hold on
    fully_fcc_idx = find(Z_eq_FCC(:,3) ~= 0);
    not_fully_fcc_idx = find(Z_eq_FCC(:,3) == 0);
    h(1) = scatter(PVD_compo(fully_fcc_idx,1), PVD_compo(fully_fcc_idx,2), 300, [0 0.4470 0.7410],'filled');
    h(2) = scatter(PVD_compo(not_fully_fcc_idx,1), PVD_compo(not_fully_fcc_idx,2), 300, [0.500 0.500 0.500],'filled');

    xlabel('position x', 'FontSize',14);
    ylabel('position y', 'FontSize',14);
    box on; axis square;
    xticks(0:10:100); xtickangle(45); yticks(0:10:100);
    set(gca,'xscale','lin', 'FontSize',14);
    set(gca,'yscale','lin', 'FontSize',14);
    hold off;

    labels = num2str((1:size(PVD_compo,1))','%d');
    text(PVD_compo(:,1), PVD_compo(:,2), labels, ...
        'horizontal','center', 'vertical','middle', ...
        'color', [0.8 0.8 0.8], ...
        'FontSize', 12, ...
    'FontWeight', 'bold')

    legend(h([1 2]), 'fully FCC', 'not fully FCC', ...
        'Interpreter', 'none', ...
        'Location', 'bestoutside', ...   
        'FontSize', 16);

    % lgd = legend;
    % lgd.Layout.Tile = 2;

    saveas(gcf, [char(group(permu_i,1)), '-',char(group(permu_i,2)), '-', char(group(permu_i,3)), ...
        '-', char(group(permu_i,4)), '-', char(group(permu_i,5)),'_Gmin-FullEquil_FullEquil', '.pdf'],'pdf');

    %% Output the palette_table

    % find the fcc spots in the Gmin diagram
    count_gm_Gmin_FCC_temp = length(find(Z_gm_Gmin_FCC(:,3)==1))/length(PVD_compo)*100;

    % find the fcc spots in the full equilibrium diagram
    count_eq_FCC_temp = length(find(Z_eq_FCC(:,3)==1))/length(PVD_compo)*100;

    % find the overlap fcc spots in both diagram
    count_gm_Gmin_eq_FCC_temp = length(find(Z_gm_Gmin_eq_FCC(:,3)==1))/length(PVD_compo)*100;


    count_fcc_fractions(permu_i,:) = [count_gm_Gmin_FCC_temp count_eq_FCC_temp count_gm_Gmin_eq_FCC_temp];
    count_ABCDE(permu_i,:)         = [char(group(permu_i,1)), '-',char(group(permu_i,2)), ...
        '-', char(group(permu_i,3)), '-', char(group(permu_i,4)), '-', char(group(permu_i,5))];

    % Create a table with the data
    palette_table = table(string(count_ABCDE), count_fcc_fractions(:,1), count_fcc_fractions(:,2), count_fcc_fractions(:,3), ...
        'VariableNames', {'ABCDE', 'Gmin_FCC_fraction', 'eq_FCC_fraction', 'Gmin_eq_FCC_fraction'});
    filename1 = 'count_fcc_fractions_byPalette.xlsx';
    % Write the table to an Excel file
    writetable(palette_table, filename1);

    %% Output the compo_table
    filename2 = 'count_fcc_fractions_byCompo.xlsx';
    writetable(compo_table, filename2);


end

 
%%









