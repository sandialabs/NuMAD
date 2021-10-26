function [model] = readInput(fstfn,bladefn,adfn,pitchAxisDomain,pitchAxisMat,LCSvalsMat)
%readInput Reads input from input data files required for analysis
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [model] = readInput(fstfn,bladefn,adfn,pitchAxisMat,LCSvalsMat)
%                    
%   This function opens and reads input files associated with analysis and
%   stores data in the model struct. The FAST file (FST) is used for 
%   reading in hub radius, tip radius, and calculating the actual blade
%   length. The blade file is used to read in geometric and structural 
%   properties of the blade. The aerodyn file (IPT) is used to read in 
%   section airfoil and chord information. This function also generates
%   boundary condition information for a cantilevered blade in the hub
%   frame.
%
%      input:
%      fstfn           = FAST filename
%      bladefn         = blade filename
%      adfn            = AeroDyn filename
%      pitchAxisMat    = filename for MATLAB .MAT containing pitch axis
%                        (reference axis) values at blade stations
%      pitchAxisDomain = spanwise location pitch axis is define at
%      LCSvalsMat      = filename for MATLAB .MAT containing lift curve
%                         slope values of blade airfoils
%
%      output:
%      model           = the model struct is output, containing all the
%                        information read in from input files that is 
%                        required for analysis


%process fast input file for data
fast=readFastMain(fstfn);
model.hubRadius=fast.TurbConf.HubRad;
model.bladeLength=fast.TurbConf.TipRad-fast.TurbConf.HubRad;

%hardwire cantilevered BC at root
model.BC.numpBC = 6;
model.BC.numsBC = 0;
model.BC.nummBC = 0;
model.BC.pBC = [1 1 0.0;
                1 2 0.0;
                1 3 0.0;
                1 4 0.0;
                1 5 0.0;
                1 6 0.0];

model.elementOrder = 1;

%load pitchAxis and LCSvals to memory
blade = readFastBlade(bladefn);
model.LCS = LCSvalsMat;
R =  blade.prop.BlFract.*model.bladeLength;

%map pitch axis from NuMAD domain to FAST blade domain
model.ptchAx = interp1(pitchAxisDomain,pitchAxisMat,R,'linear');

%Read Section Properties file
[model] = readSectionPropFile(bladefn,adfn,model.bladeLength,model);

end
