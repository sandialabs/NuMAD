% set up paths 
clc
% add path for NuMAD toolbox
NuMAD_path = 'C:\DesignCodes\NuMAD\source\';
global ANSYS_Path
ANSYS_Path = 'C:\Program Files\ANSYS Inc\v201\ansys\bin\winx64\ANSYS201.exe';
addpath(genpath(NuMAD_path))

disp('Design Code path setup script complete.  (including NuMAD, PRENUMAD, and rotor design tools)')