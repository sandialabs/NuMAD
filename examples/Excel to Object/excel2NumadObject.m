%% Example how to construct blade from excel sheet and build .yaml

close all
clear all
clc

addNumadPaths

%% Design folder


designFile = 'Excel2ObjectExample.xlsx';

blade = xlsBlade(designFile);


%% update the blade object
blade.updateBlade

%% generate a blade input file using PreComp
blade.generateBeamModel

%% generate a NuMAD file
NuMADfile = 'Excel2ObjectExample.nmd';
BladeDef_to_NuMADfile(blade,NuMADfile,'MatDBsi.txt');
numad(NuMADfile)

%% generate an ANSYS mesh
% blade.mesh = 0.1; % mesh size for ANSYS [m?]
% layupDesign_ANSYSmesh(blade,NuMADfile);