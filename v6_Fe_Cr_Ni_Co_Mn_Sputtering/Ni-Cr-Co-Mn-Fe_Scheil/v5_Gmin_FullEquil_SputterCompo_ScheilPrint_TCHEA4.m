% 20200407
% Calculation of Fe-C-Mn-Al alloy mn content in cementite at 575C

clc,clearvars,close all

% addpath('H:\Matlab Toolbox HEA\FunctionLib');
% % % add_ternary_paths

%% Define the contour
T = importdata('SputteringCompoMapNormalised.dat');
tk_delta = 298.15;
% tk      = (900+273.15) : 50 : (1300+273.15);

% outputname = ['KW130_FeCrNi-MoTi_', num2str(tk_delta),'K'];

group = [ ...
    {'Ni'}, {'Cr'}, {'Co'}, {'Mn'}, {'Fe'}; ...
    ];

%% Define contour matrix and start calculation

% for i = 1: size(group,1)
i =1;

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
        num2str(x3),'Co-', num2str(x4),'Mn-',num2str(x5),'Fe'];

    old = ".";
    new = "_";
    scheil_fname = replace(scheil_fname,old,new);

    disp(scheil_fname);

    cd('NiCrCoMnFe_Scheil_TCHEA4');
    status = copyfile('Scheil_NiCrCoMnFe_TCHEA4_template.TCM', [scheil_fname '.TCM']);

    fid = fopen([scheil_fname '.TCM']); % this file has to already be created.
    n = 1;
    textline = {};

    while( ~feof(fid) ) % This just runs until the end of the file is reached.
        textline(n,1) = {fgetl(fid)}; % This will read in every line in the file as well.

        if ( n == 8 )
            textline(n,1) = {['Ni ', num2str(x1)]};
        end

        if ( n == 9 )
            textline(n,1) = {['Cr ', num2str(x2)]};
        end

        if ( n == 10 )
            textline(n,1) = {['Co ', num2str(x3)]};
        end

        if ( n == 11 )
            textline(n,1) = {['Mn ', num2str(x4)]};
        end

        if ( n == 32 )
            % temp = cell2mat(textline(n));
            textline(n,1) = {['make-exp-data file OUTPUT_', scheil_fname, '_FCC_NS_T.exp Y']};
        end

        n = n + 1;
    end
    fclose(fid);

    fid = fopen([scheil_fname '.TCM'], 'w'); % rewriting the file
    for n = 1:length(textline)
        fprintf(fid,'%s\n', cell2mat(textline(n)));
    end

    cd ..

    % end of modify .TCM

end


