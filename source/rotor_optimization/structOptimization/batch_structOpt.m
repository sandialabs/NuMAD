% run a parallel simulation in batch mode

% set up paths
try
    snl_thor
catch
    % the paths have already been initialized
end

% rotor design model folder location
cd \\thor-storage\scratch\blennis\RotorDesign\SNL3MW\SNL3p0-148-mk0p2\NuMAD

% parallel preferences
profile_name = 'thor';
num_workers = 64;

% controller model files need to be added to the batch simulation
run('../runIEC_ipt.m')
% % directory_files = dir(params.simulinkModelFolder);
% % filenames = {directory_files(3:end).name};
% % foldernames = {directory_files(3:end).folder};
% % controller_files = cellfun('strcat(foldernames,filenames)',filenames,foldernames)
files_other_directory = {params.simulinkModelFolder};


% add in files for every type of simulation
% EXAMPLE: files_other_directory = [files_other_directory,'other_file_1.mdl','other_file_2.m'];
% files_other_directory = [files_other_directory,'\\thor-storage\scratch\blennis\RotorDesign\SNL3MW\delete--SNL3p0-148-mk0p2\'];

% create the job and run the script
job = batch('initOpt_mass_snl3p0', 'Profile',profile_name, 'Pool',num_workers-1,...
    ...'AttachedFiles',files_other_directory, ...
    'CaptureDiary',true);


