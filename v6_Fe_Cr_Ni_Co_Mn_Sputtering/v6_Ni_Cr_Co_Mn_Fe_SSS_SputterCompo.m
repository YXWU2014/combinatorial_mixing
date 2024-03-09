% 20200407

clc,clearvars,close all
addpath('H:\Matlab Toolbox HEA\FunctionLib');

%% Define the input
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

compo_table =  table([], [], [], [], [], ...
    [], [], [], [], ...
    'VariableNames', {ele_A, ele_B, ele_C, ele_D, ele_E, ...
    'ShearModulus_eff', 'PoissonsRatio_eff', 'delta_prime_misfit', 'sigma_SSS'});

% count_fcc_fractions = nan(size(group,1),3);
% count_ABCDE         = char(nan(size(group,1), length([char(group(1,1)), '-',char(group(1,2)), ...
%     '-', char(group(1,3)), '-', char(group(1,4)), '-', char(group(1,5))])));

%%

for permu_i = 1: size(group,1)

    %     % Skip the calculated
    %     if isfile([char(group(permu_i,1)), '-',char(group(permu_i,2)),'-', char(group(permu_i,3)),...
    %             '-', char(group(permu_i,4)), '-', char(group(permu_i,5)), '_Sputter compo contour_', num2str(tk_Gmin),'K_','.jpeg'])
    %         continue;
    %     end

    % permu_i = 1

    %% Initiate TC system
    tc_init_root;
    tc_open_database('tchea4');
    tc_check_error;

    % Select elements and retrieve data

    ele_A_permu_i = char(group(permu_i,1));
    ele_B_permu_i = char(group(permu_i,2));
    ele_C_permu_i = char(group(permu_i,3));
    ele_D_permu_i = char(group(permu_i,4));
    ele_E_permu_i = char(group(permu_i,5));

    sel = [ele_A_permu_i, ' ', ele_B_permu_i, ' ',ele_C_permu_i, ' ', ...
        ele_D_permu_i, ' ', ele_E_permu_i];
    tc_element_select(sel);
    tc_get_data;

    % Set conditions
    tc_set_condition('n',1);
    tc_set_condition('p',101325);


    %% tabulated data from pymatgen and TC

    %         ele_A = 'V';
    %         ele_B = 'Co';
    %         ele_C = 'Ni';
    %         ele_D = 'Fe';
    %         ele_E = 'Cr';

    % from pymatgen
    ele_all           = {'Fe', 'Cr', 'Ni', 'Co', 'V', 'Mn', 'Mo', ...
        'Cu', 'Nb', 'W', 'Ti', 'Al', 'Si', 'Ta'};
    molarV_matgen_ele = {7.09  7.23  6.59  6.67  8.32  7.35  9.38 ...
        7.11  10.83  9.47 10.64 10.00 12.06 10.85};
    YM_matgen_ele     = {211.0 279.0 200.0 209.0 128.0 198.0 329.0 ...
        130.0 105.0 411.0 116.0  70.0  47.0 186.0};
    BM_matgen_ele     = {170.0 160.0 180.0 180.0 160.0 120.0 230.0 ...
        140.0 170.0 310.0 110.0  76.0 100.0 200.0};

    % from TC
    AtomicV_FCC_ele = {11.375143, 12.441593, 10.932037, 11.122192, 14.685831, 12.665726, 15.924223, ...
        11.794933, 17.893725, 16.509974, 17.647633, 16.570502, 14.431877, 18.140344};

    % Find the indices of the elements in ele_matgen
    %     ele_indices = find(contains(ele_all, {ele_A_permu_i, ele_B_permu_i, ele_C_permu_i, ele_D_permu_i, ele_E_permu_i}));
    ele_indices = [find(contains(ele_all, {ele_A_permu_i})), ...
        find(contains(ele_all, {ele_B_permu_i})), ...
        find(contains(ele_all, {ele_C_permu_i})), ...
        find(contains(ele_all, {ele_D_permu_i})), ...
        find(contains(ele_all, {ele_E_permu_i}))];


    % Retrieve the corresponding values for each element
    molarV_matgen = molarV_matgen_ele(ele_indices);
    YM_matgen = YM_matgen_ele(ele_indices);
    BM_matgen = BM_matgen_ele(ele_indices);
    AtomicV_FCC = AtomicV_FCC_ele(ele_indices);

    % Assign each value to a variable with the corresponding element name
    [molarV_matgen_A, molarV_matgen_B, molarV_matgen_C, molarV_matgen_D, molarV_matgen_E] = deal(molarV_matgen{:});
    [YM_matgen_A, YM_matgen_B, YM_matgen_C, YM_matgen_D, YM_matgen_E]                     = deal(YM_matgen{:});
    [BM_matgen_A, BM_matgen_B, BM_matgen_C, BM_matgen_D, BM_matgen_E]                     = deal(BM_matgen{:});
    [AtomicV_FCC_A, AtomicV_FCC_B, AtomicV_FCC_C, AtomicV_FCC_D, AtomicV_FCC_E]           = deal(AtomicV_FCC{:});
    molarV_ABCDE      = [molarV_matgen_A, molarV_matgen_B, molarV_matgen_C, molarV_matgen_D, molarV_matgen_E];
    YM_ABCDE          = [YM_matgen_A, YM_matgen_B, YM_matgen_C, YM_matgen_D, YM_matgen_E];
    BM_ABCDE          = [BM_matgen_A, BM_matgen_B, BM_matgen_C, BM_matgen_D, BM_matgen_E];

    YM_eff = NaN(length(PVD_compo),1);
    BM_eff = NaN(length(PVD_compo),1);

    for j = 1: length(PVD_compo)
        PVD_compo_ABCDE_j = PVD_compo(j, 3:7);
        YM_eff(j)         = sum(PVD_compo_ABCDE_j .* molarV_ABCDE .* YM_ABCDE) / sum(PVD_compo_ABCDE_j.* molarV_ABCDE);
        BM_eff(j)         = sum(PVD_compo_ABCDE_j .* molarV_ABCDE .* BM_ABCDE) / sum(PVD_compo_ABCDE_j.* molarV_ABCDE);
    end

    nu_eff = (1-YM_eff./3./BM_eff)/2; % Poisson's ratio by isotropic elasticity
    SM_eff = YM_eff/2./(1+nu_eff);   % effective shear modulus

    %% --- calculate delta prime misfit volume

    tk_delta = 298.15;

    % Initialize the output arrays
    n_compositions = length(PVD_compo);

    Z_LP_FCC_L12         = [PVD_compo(:,[1 2]), NaN(length(PVD_compo),1)];
    Z_VA_FCC_L12         = [PVD_compo(:,[1 2]), NaN(length(PVD_compo),1)];
    Z_b_FCC_L12          = [PVD_compo(:,[1 2]), NaN(length(PVD_compo),1)];
    Z_delta_prime_misfit = [PVD_compo(:,[1 2]), NaN(length(PVD_compo),1)];

    % Loop over each composition
    for j = 1:n_compositions

        % Set the composition and temperature conditions
        tc_set_condition(['x(', char(group(permu_i,1)), ')'], PVD_compo(j,3));
        tc_set_condition(['x(', char(group(permu_i,2)), ')'], PVD_compo(j,4));
        tc_set_condition(['x(', char(group(permu_i,3)), ')'], PVD_compo(j,5));
        tc_set_condition(['x(', char(group(permu_i,4)), ')'], PVD_compo(j,6));

        tc_set_condition('T', tk_delta);
        tc_set_phase_status('*', 'suspended', 0); %@/0/or/1/Start value, number of mole formula units
        tc_set_phase_status('FCC_L12', 'entered', 1);
        tc_set_minimization('off')

        tc_compute_equilibrium;

        [a,b]=tc_error;
        if a~=0
            s=sprintf('ERROR %d',a);
            disp(s);
            disp(b);
            gm_FCC_L12_temp = nan;
            vm_FCC_L12_temp = nan;
        else
            vm_FCC_L12_temp = tc_get_value('vm(FCC_L12)');
            LP_FCC_L12_temp = ((4*vm_FCC_L12_temp/6.02214179E23)^(1/3))*1E10; % in A
            VA_FCC_L12_temp = (vm_FCC_L12_temp/6.02214179E23)*1E30; % in A

            % average elemental misfit volume (A^3)
            delta_V_A_temp = abs(AtomicV_FCC_A - VA_FCC_L12_temp);
            delta_V_B_temp = abs(AtomicV_FCC_B - VA_FCC_L12_temp);
            delta_V_C_temp = abs(AtomicV_FCC_C - VA_FCC_L12_temp);
            delta_V_D_temp = abs(AtomicV_FCC_D - VA_FCC_L12_temp);
            delta_V_E_temp = abs(AtomicV_FCC_E - VA_FCC_L12_temp);

            % average Burgers vector for CCA (A)
            b_FCC_L12_temp = ((4*VA_FCC_L12_temp)^(1/3))/(2^(1/2));

            % patch composition
            c_A_temp = PVD_compo(j,3);
            c_B_temp = PVD_compo(j,4);
            c_C_temp = PVD_compo(j,5);
            c_D_temp = PVD_compo(j,6);
            c_E_temp = 1- c_A_temp - c_B_temp - c_C_temp - c_D_temp;

            % definition of misfit parameter delta_prime in Nohring2019
            delta_prime_misfit_temp = sqrt(c_A_temp * delta_V_A_temp^2 + c_B_temp * delta_V_B_temp^2 + ...
                c_C_temp * delta_V_C_temp^2 + c_D_temp * delta_V_D_temp^2 + ...
                c_E_temp * delta_V_E_temp^2)/b_FCC_L12_temp^3;

            % save into the contour matrix
            Z_LP_FCC_L12(j,3)         = LP_FCC_L12_temp;
            Z_VA_FCC_L12(j,3)         = VA_FCC_L12_temp;
            Z_b_FCC_L12(j,3)          = b_FCC_L12_temp;
            Z_delta_prime_misfit(j,3) = delta_prime_misfit_temp;

        end
        tc_reset_error;
    end

    %% calculate the solid solution contribution

    kb    = 1.380649e-23;        % Boltzmann constant J*K−1
    M     = 3.06;                % Taylor factor
    alpha = 0.123;               % The line tension parameter obtained from edge dislocation line tension in the EAM FeNiCr effective matrix

    mu_eff = SM_eff * 1E9; % effective shear modulus in [Pa]
    % nu_eff  % effective Poisson's ratio

    % Peierls Stress at 0K
    tao_0   = 0.051*0.35*alpha^(-1/3) * mu_eff .* ((1+nu_eff)./(1-nu_eff)).^(4/3) .* (Z_delta_prime_misfit(:,3).^2).^(2/3)/1E6; % in MPa
    sigma_0 = tao_0 * M; % in MPa

    % Energy barrier
    delta_Eb = 0.274*5.7*alpha^(1/3) * mu_eff .* (Z_b_FCC_L12(:,3).^3)*1e-30 .* ((1+nu_eff)./(1-nu_eff)).^(2/3) .* (Z_delta_prime_misfit(:,3).^2).^(1/3); % in J
    % disp([num2str(delta_Eb*6.241509E18) ' ev']);

    epsilon0 = 1E4;  % reference strain rate s-1
    epsilon  = 1E-3; % strain rate s-1

    sigma_T = sigma_0 .* (1 - (kb*tk_delta./delta_Eb * log(epsilon0/epsilon)).^(2/3)); % in MPa

    %% extract the SSS parameters

    compo_data_j = table(PVD_compo(:, 3), PVD_compo(:, 4), PVD_compo(:, 5), PVD_compo(:, 6), PVD_compo(:, 7),...
        SM_eff, nu_eff, Z_delta_prime_misfit(:, 3), sigma_T, ...
        'VariableNames', {char(group(permu_i,1)), char(group(permu_i,2)), char(group(permu_i,3)), char(group(permu_i,4)), char(group(permu_i,5)), ...
        'ShearModulus_eff', 'PoissonsRatio_eff', 'delta_prime_misfit', 'sigma_SSS'});

    compo_data_j = movevars(compo_data_j, ele_B, 'After', ele_A);
    compo_data_j = movevars(compo_data_j, ele_C, 'After', ele_B);
    compo_data_j = movevars(compo_data_j, ele_D, 'After', ele_C);
    compo_data_j = movevars(compo_data_j, ele_E, 'After', ele_D);

    compo_table = vertcat(compo_table, compo_data_j);

    %% plot - sss
    close all;
    figure(1);
    cmap = linspecer;
    hold on
    set(gcf,'units','points','position',[0,0,600,250]);
    subplot(1, 1, 1)
    scatter(PVD_compo(:,1), PVD_compo(:,2), 300, sigma_T,"filled")
    scatter(PVD_compo(:,1), PVD_compo(:,2), 300, [0.5 0.5 0.5],'LineWidth', 1)
    colormap(gca, cmap)

    c= colorbar;
    % caxis([100 400]);
    % title('Misfit Parameter($\delta^{\prime}$)','Interpreter', 'latex');
    xlabel('position x','FontSize',14);
    ylabel('position y','FontSize',14);
    c.FontSize = 16;
    c.Label.String = 'Solid solution strengthening (MPa)';
    c.Label.FontSize = 16;

    box on; axis square; xticks(0:10:100); yticks(0:10:100);
    set(gca,'xscale','lin', 'FontSize',14); xtickangle(45);
    set(gca,'yscale','lin', 'FontSize',14);

    labels = num2str((1:size(PVD_compo,1))','%d');
    text(PVD_compo(:,1), PVD_compo(:,2), labels, ...
        'horizontal','center', 'vertical','middle', ...
        'color', [0.5 0.5 0.5], ...
        'FontSize', 12, ...
        'FontWeight', 'bold')

    saveas(gcf, [char(group(permu_i,1)), '-',char(group(permu_i,2)), '-', char(group(permu_i,3)), ...
        '-', char(group(permu_i,4)), '-', char(group(permu_i,5)),'_SolidSolutionStrengthening_', num2str(tk_delta),'K', '.pdf'],'pdf');

    %% Output the compo_table
    filename2 = 'SSS_byCompo.xlsx';
    writetable(compo_table, filename2);

end






%%




