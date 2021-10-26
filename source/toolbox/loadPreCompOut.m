function out = loadPreCompOut(varargin)
% LOADPRECOMPOUT  Load PreComp output data, *_gen files only
%                          Under Construction
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% Requires: txt2mat.m (MatlabCentral File ID: #18430)
%

if nargin == 0
    % filter files with .out extension
    arg = '*.out_gen';
else
    % pass all arguments
    arg = varargin{:};
end
[A,ffn,nh,SR,hl,fpos] = txt2mat(arg);
% return if file selection cancelled
if isempty(ffn); return; end;
% parse header lines, assuming labels and units are at end of header
head = regexp(hl,'[^\n]+','match');
out.list = strtrim(regexp(head{end-1},'[^\s]+','match'));
out.units = strtrim(regexp(head{end},'[^\s]+','match'));
out.data = A;
