%% Example how to construct blade from excel sheet and build .yaml

close all
clear all
clc

addNumadPaths

%% Design folder

designFolder = 'C:\Users\rclarke\Desktop\Projects\BAR\NuMAD_Test\NuMAD_Excel';
designFile = 'NREL5MW_SNL61p5.xlsx';

blade = xlsBlade(fullfile(designFolder,designFile));


%% update the blade object
blade.updateBlade


%% generate a NuMAD file
NuMADfile = 'NREL5MW_SNL61p5.nmd';
BladeDef_to_NuMADfile(blade,NuMADfile,'MatDBsi.txt');
numad(NuMADfile)

%% generate an ANSYS mesh
blade.mesh = 0.1; % mesh size for ANSYS [m?]
layupDesign_ANSYSmesh(blade,'NREL5MW_SNL61p5.nmd');