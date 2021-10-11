function [] = writeDLCoutput(output, filename, varargin)
%writeDLCoutput  write the DLC output structure to an Excel file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************

names = fieldnames(output);
xlsData = [];
for ii = 1:length(names)
    %     header = {'** DLC Name **' ['** ' names{ii} ' **'] '** Occurring at time (sec) **'...
    %         '** on Channel **' '** in File **'};
    if isempty(varargin)
        header = {'DLC Name' names{ii} 'Occurring at time (sec)' 'on Channel' 'in File'};
        xlsData = [xlsData; header; table2cell(output.(names{ii}))]; %#ok<AGROW>
    elseif strcmp(varargin, '1')
        header = {'DLC Name' names{ii} 'Occurring at time (sec)' 'on Channel' 'in File' 'Variable Name'};
        xlsData = [xlsData; header; table2cell(output.(names{ii}))];   
    elseif strcmp(varargin, '2')
        header = {'Theta' 'DLC Name' names{ii} 'Occurring at time (sec)' 'on Channel' 'in File'};
        xlsData = [xlsData; header; [output.(names{ii}).Properties.RowNames table2cell(output.(names{ii}))]]; 
    elseif strcmp(varargin, '3')
        header = {names{ii} 'Fx' 'Fy' 'Fz' 'Mx' 'My' 'Mz'};
        xlsData = [xlsData; header; [output.(names{ii}).Properties.RowNames table2cell(output.(names{ii}))]]; 
    end    
end

xlswrite(filename, xlsData,1);