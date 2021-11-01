function settings = readNuMADsettings(filename)
%READNUMADSETTINGS  Read the NuMAD settings file 
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   settings = readNuMADsettings(FILENAME)
%     Read the NuMAD settings file FILENAME (typically 'settings.txt' 
%     located in %APPDATA%\NuMAD\) and return data structure.
%
%   See also writeNuMADsettings

%===== CREDITS & CHANGELOG ================================================
%Developed by Wind & Water Power Technologies, Sandia National Laboratories
%2011.02.16  JCB: first draft, based on readNuMADinput()

% ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% don't read settings file when in parallel
workID = getCurrentWorker;
global ansysPath
global bmodesPath
global precompPath
if isempty(workID) % not running on parallel node/worker - read settings
    
    
    %filename = 'NuMAD_settings.txt';
    
    % Open the file and read the entire contents
    fid = fopen(filename);
    if (fid == -1)
        error('Could not open file "%s"',filename);
    end
    filecontents = fread(fid,inf,'uint8=>char')';
    fclose(fid);
    
    
    % %%%%%%%%%%%%  PROGRAMMER'S NOTE %%%%%%%%%%%%%
    % The following regular expression pattern matches any number of characters
    % found between the opening and closing "reference" tags
    %pattern = '<reference>(.*)</reference>';
    %t = regexp(filecontents, pattern, 'tokens');
    % t is a cell containing a cell array
    % try
    %     % jcb: is there a better way to extract the contents of t?
    %     af.reference = cell2mat(t{1});
    % catch me
    %     af.reference = '';
    % end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pattern = '<ansys_path>(.*)</ansys_path>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {''};
    end
    settings.ansys_path = cell2mat(t{1});
    
    pattern = '<ansys_product>(.*)</ansys_product>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {''};
    end
    settings.ansys_product = cell2mat(t{1});
    
    pattern = '<bmodes_path>(.*)</bmodes_path>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {''};
    end
    settings.bmodes_path = cell2mat(t{1});
    
    pattern = '<precomp_path>(.*)</precomp_path>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {''};
    end
    settings.precomp_path = cell2mat(t{1});
    
    pattern = '<gui_layout monitors="([^"]*)">';
    t = regexp(filecontents, pattern, 'tokens');
    settings.monitors = t{1}{1};
    
    pattern = '<xy_main>(.*)</xy_main>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f %f %f');
    settings.xy_main = cell2mat(values);
    
    pattern = '<xy_materials>(.*)</xy_materials>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f');
    settings.xy_materials = [values{1} values{2}];
    
    pattern = '<xy_isotropic>(.*)</xy_isotropic>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f');
    settings.xy_isotropic = [values{1} values{2}];
    
    pattern = '<xy_orthotropic>(.*)</xy_orthotropic>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f');
    settings.xy_orthotropic = [values{1} values{2}];
    
    pattern = '<xy_composite>(.*)</xy_composite>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f');
    settings.xy_composite = [values{1} values{2}];
    
    pattern = '<xy_activelist>(.*)</xy_activelist>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f');
    settings.xy_activelist = [values{1} values{2}];
    
    pattern = '<xy_ansys>(.*)</xy_ansys>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f');
    settings.xy_ansys = [values{1} values{2}];
    
    pattern = '<xy_genline>(.*)</xy_genline>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f');
    settings.xy_genline = [values{1} values{2}];
    
    pattern = '<xy_plot3d>(.*)</xy_plot3d>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f');
    settings.xy_plot3d = [values{1} values{2}];
    
    pattern = '<xy_flutterinput>(.*)</xy_flutterinput>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f');
    settings.xy_flutterinput = [values{1} values{2}];
    
    pattern = '<xy_pathconfig>(.*)</xy_pathconfig>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50'};
    end
    values = textscan(cell2mat(t{1}),'%f %f');
    settings.xy_pathconfig = [values{1} values{2}];
    
    pattern = '<xy_bpesegments>(.*)</xy_bpesegments>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {'20 50 800 300'};
    end
    values = textscan(cell2mat(t{1}),'%f %f %f %f');
    settings.xy_bpesegments = [values{1} values{2} values{3} values{4}];
    
    pattern = '<job_path>(.*)</job_path>';
    t = regexp(filecontents, pattern, 'tokens');
    if isempty(t)
        t{1} = {''};
    end
    settings.job_path = cell2mat(t{1});
    
    % DEBUGGING
    %assignin('base','settings',settings);
    
else % running on parallel worker -- don't read settings file
    % define settings file for parallel operation    
    disp('Parallel Operation: do not read NuMAD settings file...')

    settings.ansys_path = ansysPath;
    settings.ansys_product = 'ANSYS';
    settings.bmodes_path = bmodesPath;
    settings.precomp_path = precompPath;
    settings.monitors = '';
    settings.xy_main = [20 50 1000 600];
    settings.xy_materials = [20 50];
    settings.xy_isotropic = [20 50];
    settings.xy_orthotropic = [20 50];
    settings.xy_composite = [20 50];
    settings.xy_activelist = [20 50];
    settings.xy_ansys = [20 50];
    settings.xy_genline = [20 50];
    settings.xy_plot3d = [20 50];
    settings.xy_pathconfig = [20 50];
    settings.xy_flutterinput = [20 50];
    settings.xy_bpesegments = [20 50 800 300];
    if ispc
        settings.job_path = getenv('USERPROFILE');
    elseif isunix
        settings.job_path = getenv('HOME');
    else
        settings.job_path = '';
    end
    settings.job_name = '';
    
end
% ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
