clear all
clc

% cd('C:\data\_BAR\test')

% Create blade in prenumad
blade = BladeDef;
blade.readYAML('BAR00.yaml')
return
changeBlade = 0;
if changeBlade
    % Change blade object dimensions
    % Example make both shear web cores 10% thicker
    blade.components(13).cp(:,2) = blade.components(13).cp(:,2).*1.1;
    blade.components(16).cp(:,2) = blade.components(16).cp(:,2).*1.1;
end

blade.updateBlade

% create NuMAD input file for NuMAD GUI
% hm = pwd;
% cd('NuMAD')
BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt');

% generate ANSYS mesh
movefile('yaml2BladeDef_BAR008_ble/airfoils')
movefile('yaml2BladeDef_BAR008_ble/BAR008_ble_blade.mat')
numad('numad.nmd','ansys',blade.mesh)

% generate a FAST beam model from PreComp and BModes
% cd(hm)
blade.establishPaths
blade.generateBeamModel

save blade blade
if changeBlade
    % Write modified blade object to yaml file
    blade.writeYAML('BAR008_mod.yaml')
end