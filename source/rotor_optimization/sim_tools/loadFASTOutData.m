function out = loadFASTOutData(varargin)
% LOADFASTOUTDATA  Read FAST OutData into structure
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% Requires: txt2mat.m (MatlabCentral File ID: #18430)
%
% Syntax:
%   out = loadFASTOutData;
%   out = loadFASTOutData('filename');
%
%   where structure 'out' has fields {'data','list','units'}
%
% Plotting:
%   plot(out.data(:,1),out.data(:,N));
%   xlabel([out.list{1},' ',out.units{1}]); 
%   ylabel([out.list{N},' ',out.units{N}]);
% 
%   See also plotFASTOutData

if nargin == 0
    % filter files with .out extension
    arg = '*.out';
else
    % pass all arguments
    arg = varargin{:};
end

[A,ffn,nh,SR,hl,fpos] = txt2mat(arg);
% return if file selection cancelled
if isempty(ffn); return; end;
% parse header lines, assuming labels and units are at end of header
head = regexp(hl,'[^\n]+','match');
out.list = strtrim(regexp(head{end-1},'[^\t]+','match'));
out.units = strtrim(regexp(head{end},'[^\t]+','match'));
out.data = A;
% ble: create a table format as well -- easier to use in MATLAB
out.tbl = array2table(A);
out.tbl.Properties.VariableNames = out.list;


