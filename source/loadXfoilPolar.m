function out = loadXfoilPolar(varargin)
% LOADXFOILPOLAR  Read XFOIL results polar (filename.txt) into structure
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% Requires: txt2mat.m (MatlabCentral File ID: #18430)
%
% Syntax:
%   out = loadXfoilPolar;
%   out = loadXfoilPolar('filename');
%
%   where structure 'out' has fields {'data','list'}
%
% Plotting:
%   plot(out.data(:,1),out.data(:,N));
%   xlabel(out.list{1}); 
%   ylabel(out.list{N});
% 

if nargin == 0
    % filter files with .out extension
    arg = '*.dat';
else
    % pass all arguments
    arg = varargin{:};
end
[A,ffn,nh,SR,hl,fpos] = txt2mat(arg);
% return if file selection cancelled
if isempty(ffn); return; end;
% parse header lines, assuming labels and units are at end of header
head = regexp(hl,'[^\n]+','match');
out.title= strtrim(head{4});
a=regexp(head{9},'[\w.]+','match');
out.ReNum= str2num(a{4})*10^str2num(a{6});
out.list = strtrim(regexp(head{11},'\w+','match'));
out.data = A;
