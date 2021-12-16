% add path for NuMAD source and optional toolbocxes
clc
global numadPath
global ansysPath
global bmodesPath
global precompPath
global fastPath
global crunchPath
global adamsPath  
global turbsimPath
global iecwindPath

numadPath = fullfile(pwd,'source');
ansysPath = 'C:\Program Files\ANSYS Inc\v201\ansys\bin\winx64\ANSYS201.exe';
precompPath = 'C:\DesignCodes\PreComp_v1.00.03\PreComp.exe';
bmodesPath = 'C:\DesignCodes\BModes_v3.00.00\BModes.exe';
fastPath = 'C:\DesignCodes\FAST_v7.02.00d\FAST.exe';
crunchPath = 'C:\DesignCodes\Crunch_v3.00.00\Crunch.exe';
turbsimPath='C:\DesignCodes\TurbSim_v1.50\TurbSim.exe';
iecwindPath='C:\DesignCodes\IECWind\IECWind.exe';

addpath(genpath(numadPath))

disp('NuMAD and Design Code path setup complete.')