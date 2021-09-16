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
%   Please see the NuMAD User's Guide for detailed usage instructions.
%

    % this function simply serves as an alias and wrapper for NuMAD_main
    NuMAD_main(varargin{:});
end