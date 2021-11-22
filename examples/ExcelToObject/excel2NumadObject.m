%% Example how to construct blade from excel sheet and build .yaml
% Verify the NuMAD paths are added prior to running this example.
% This can be done by executing the ``addNumadPaths`` script, shown below.

<<<<<<< HEAD
close all
clc

addNumadPaths
=======
run('../../addNumadPaths')
>>>>>>> 61b38be64c420aef0ae10bac117569cbf17d0204

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