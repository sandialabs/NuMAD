function generateBLASTFreqDampPlots(matfile)
%generatBLASTFreqDampPlots  generates freq/damp plots from saved data
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   generateBLASTFreqDampPlots(matfile)
%                    
%   This functions loads frequency, damping, and operating condition
%   information from a saved .mat file and generates freqeuncy vs. op.
%   condition and damp vs. op condition plots
%
%      input:
%      matfile      = name of .mat file containing freq/damp information
 
%      output:
%      none

load(matfile);
OmegaArray = omegaPlot;
plotFreqDampCurves(OmegaArray,convergedDamp,convergedFreq,analysisType);

end

