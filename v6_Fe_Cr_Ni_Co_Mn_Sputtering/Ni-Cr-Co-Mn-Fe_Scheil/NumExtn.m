
% Yuxiang WU 20170330
% Extract the number data from DICTRA OUTPUT (i.e. remove all the text)

function Data = NumExtn(fname)
Data = [];

for j=1:length(fname) % OUTPUTS to be analysed
    
    Data_temp = [];
    finput = fopen(char(fname(j)),'rt');
    
    if finput == 0
        error('CANNOT read the OUTPUT file')
    end
    
    % avoid  Err "Error using str2num (line 31) Requires character vector or array input."
    %l = fgetl(finput);
    %while ~ischar(l) || ~ismatrix(l)
    %    pause(0.1);
    %    l = fgetl(finput);
    %end
    
    % Read line-by-line
    while ~feof(finput)
        
        l = fgetl(finput);
        
        % this is to remove the 'M' in the first line of output.exp
        l(l=='M')='';
        
        % tline is a character vector unless the line contains only the end-of-file marker. 
        % In this case, tline is the numeric value -1.
        if l == -1
            continue
        end
        
        temp = str2num(l); % all the lines with string will be empty
        
        % Determine whether array is empty
        if isempty(temp)
            continue;
        end
        Data_temp = [Data_temp; temp]; % rewrite the data row-by-row
    end
    
    fclose(finput);
    
    Data = [Data Data_temp]; % Two-colume Data_temp will be stacked to the RHS
end



