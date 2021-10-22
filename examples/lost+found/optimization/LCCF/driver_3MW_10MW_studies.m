clear all; close all; clc

snl2

%% run 10 MW optimization studies
bladeFile='..\bladeObject_IEA10p0-198-mk0p4.mat';

% s2-baselineCF optimization
designNuMADfolder = 'C:\data\IEA10MW\IEA10p0-198-mk0p4-s2-baselineCF\NuMAD';
initOpt_mass_iea10p0_allCases_noSparCapWidth(designNuMADfolder, bladeFile)

% s1-fiberglass optimization
designNuMADfolder = 'C:\data\IEA10MW\IEA10p0-198-mk0p4-s1-fiberglass\NuMAD';
initOpt_mass_iea10p0_allCases_noSparCapWidth(designNuMADfolder, bladeFile)


% % s3-heavytextileCF optimization
% designNuMADfolder = 'C:\data\IEA10MW\IEA10p0-198-mk0p4-s3-heavyCF\NuMAD';
% initOpt_mass_iea10p0_allCases_noSparCapWidth(designNuMADfolder, bladeFile)


%% run 3 MW optimization studies
bladeFile='..\bladeObject_SNL3p0-148-mk0p2.mat';

% % s1-fiberglass optimization
% designNuMADfolder = 'C:\data\SNL3MW\SNL3p0-148-mk0p2-s1-fiberglass\NuMAD';
% initOpt_mass_snl3p0_allCases_noSparCapWidth(designNuMADfolder, bladeFile)
% 
% % s2-baselineCF optimization
% designNuMADfolder = 'C:\data\SNL3MW\SNL3p0-148-mk0p2-s2-baselineCF\NuMAD';
% initOpt_mass_snl3p0_allCases_noSparCapWidth(designNuMADfolder, bladeFile)
% 
% s3-heavytextileCF optimization
% designNuMADfolder = 'C:\data\SNL3MW\SNL3p0-148-mk0p2-s3-heavyCF\NuMAD';
% initOpt_mass_snl3p0_allCases_noSparCapWidth(designNuMADfolder, bladeFile)
