% run a parallel simulation in batch mode

% set up paths
try
    snl_thor
catch
    % the paths have already been initialized
end

% rotor design model folder location
% % cd \\thor-storage\SCRATCH\blennis\RotorDesign\NRT_537kg_44rpm
cd \\thor-storage\SCRATCH\blennis\RotorDesign\SWiFT_V27_v2.00.00a

% parallel preferences
profile_name = 'thor';
num_workers = 256;

% controller model files need to be added to the batch simulation
runIEC_ipt
% % directory_files = dir(params.simulinkModelFolder);
% % filenames = {directory_files(3:end).name};
% % foldernames = {directory_files(3:end).folder};
% % controller_files = cellfun('strcat(foldernames,filenames)',filenames,foldernames)
files_other_directory = {params.simulinkModelFolder};


% add in files for every type of simulation
% EXAMPLE: files_other_directory = [files_other_directory,'other_file_1.mdl','other_file_2.m'];

% create the job and run the script
job = batch('runIEC(''all'',1,0)', 'Profile',profile_name, 'Pool',num_workers-1,...
    'AttachedFiles',files_other_directory, 'CaptureDiary',true);


