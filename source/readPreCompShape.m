function output=readPreCompShape(filename)
% READPRECOMPSHAPE  Read PreComp shape file data for inspection
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% output=readPreCompShape(filename)
% Works with PreComp v1.00.03
% Requires: txt2mat.m (MatlabCentral File ID: #18430)
%
%   readPreCompShape() will open a file selector window.  Choose a PreComp
%               shape *.inp file.
%   For shape inspection only.
%   Not for use in conjunction with writePreCompShape()
%

if nargin == 0
    % open file selector and load data
    [A,ffn,nh,SR,hl,fpos] = txt2mat('*.inp');
    % return if file selection cancelled
    if isempty(ffn); return; end;
    % extract the filename
    fn =  strtrim(regexp(ffn,'[^\\]+','match'))
else
    [A,ffn,nh,SR,hl,fpos] = txt2mat(filename);
    % extract the filename
    fn =  strtrim(regexp(ffn,'[^\\]+','match'))
 end

output=A;

end

