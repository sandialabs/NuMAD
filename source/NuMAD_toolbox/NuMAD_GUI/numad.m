function numad(varargin)
%NUMAD  Launches NuMAD - a design tool for wind turbine blades
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NUMAD, called without any argument, opens the main user interface
%   
%   NUMAD(FILENAME) opens a project for editing
%      FILENAME  = filename of a NuMAD project file
%
%   NUMAD(FILENAME,batchRunType,batchArgument) opens a project for a batch run
%   with specified batchRunType:
%      'generate' - Generate ANSYS output according to saved settings
%      'shell7' - Generate only the shell7 output
%      'ansys' - Generate both the shell7 and ansysdb
%      batchArgument is the mesh element size (optional input)
%
%      'precomp' - Run PreComp-to-FASTBlade
%      batchArgument specifies which mode results are [flp1,flp2,edg1] and
%                    is a required input; typically the order is [1,3,2]
%   
%   Please see the NuMAD User's Guide for detailed usage instructions.
%

    % this function simply serves as an alias and wrapper for NuMAD_main
    NuMAD_main(varargin{:});
end