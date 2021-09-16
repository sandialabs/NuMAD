function FlutterToolInput(varargin)
%FlutterToolInput  User interface to flutter tool
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   FlutterToolInput(cbo), in standard usage, is used as a menu callback 
%   within NuMAD_main which launches the flutter tool user interface. 
%   Usage: uimenu('label','Flutter Tool','callback',@FlutterToolInput)
%
%   FlutterToolInput(FILENAME) opens flutter tool interface outside NuMAD
%      FILENAME  = filename of a NuMAD project file

%===== INITIALIZATION =====================================================

 uiwait(msgbox('The Sandia flutter tool capability has been deactivated in the current release of NuMAD.  Comparisons of flutter predictions for typical utility scale rotors (up to 60m) and very large rotors (>100m) raise important and currently unresolved questions about the approach that is currently implemented.  Sandia will release a new version of this module when the computational issues have been resolved.','Flutter analysis capability deactivated','modal'));
 
 