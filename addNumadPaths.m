% add NuMAD source and ANSYS *.exe to MATLAB path
clc

% add path for NuMAD toolbox
NuMAD_path = 'C:\DesignCodes\NuMAD\source\';
addpath(genpath(NuMAD_path))


% add path for ANSYS
global ANSYS_Path
ANSYS_Path = 'C:\Program Files\ANSYS Inc\v201\ansys\bin\winx64\ANSYS201.exe';

disp('NuMAD source added to MATLAB path, including: objects, optimization, and toolbox')