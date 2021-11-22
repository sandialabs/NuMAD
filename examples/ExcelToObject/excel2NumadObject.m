%% Example how to construct blade from excel sheet and build .yaml

close all
clc

addNumadPaths

%% Load blade object from excel sheet
designFile = 'Excel2ObjectExample.xlsx';
blade = xlsBlade(designFile);

%% update the blade object
blade.updateBlade

%% generate a blade input file using PreComp
% blade.generateBeamModel   %If PreComp is installed you can run this command to generate beam model parameters

%% generate a NuMAD file
NuMADfile = 'Excel2ObjectExample.nmd';
BladeDef_to_NuMADfile(blade,NuMADfile,'MatDBsi.txt');
numad(NuMADfile)

%% generate an ANSYS mesh
blade.mesh = 0.1; % mesh size for ANSYS 
layupDesign_ANSYSmesh(blade,NuMADfile);