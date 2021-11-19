% add path for NuMAD source and optional toolbocxes
clc
global numadPath
global ansysPath
global bmodesPath
global precompPath

numadPath = fullfile(pwd,'source');
ansysPath = 'C:\Program Files\ANSYS Inc\v181\ansys\bin\winx64\ANSYS181.exe';
precompPath = 'C:\DesignCodes\PreComp_v1.00.03\PreComp.exe';
bmodesPath = 'C:\DesignCodes\BModes_v3.00.00\BModes.exe';

addpath(genpath(numadPath))

disp('NuMAD and Design Code path setup complete.')