function snl
% User needs to enter path to their own installation of the Toolbox below
% before running this script
SNL_NuMAD_Toolbox_path = 'C:\Users\brresor\Documents\SNL_NuMAD_ToolBox';

addpath(SNL_NuMAD_Toolbox_path)
addpath(fullfile(SNL_NuMAD_Toolbox_path,'BPE'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'Docs'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'Examples'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'Examples\Materials'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'Examples\NuMAD2BPE2FASTBlade'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'Examples\NuMAD2PreComp2FASTBlade'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'Examples\NuMADforExcel'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'Examples\WTPerf'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'MatlProps'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'NuMAD'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'NuMAD\NuMAD2BPE2FASTBlade'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'NuMAD\NuMAD2PreComp2FASTBlade'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'NuMAD\airfoils'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'NuMAD\macros'))
addpath(fullfile(SNL_NuMAD_Toolbox_path,'Sub'))
p1=genpath(fullfile(SNL_NuMAD_Toolbox_path,'Docs'));
p2=genpath(fullfile(SNL_NuMAD_Toolbox_path,'Examples'));
rmpath(p1,p2);
disp('SNL NuMAD Toolbox v13.04.03 path setup script complete.')
end