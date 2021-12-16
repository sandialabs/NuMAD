function plotFASTOutData(arg)
% PLOTFASTOUTDATA  Interactive plot of FAST OutData
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% Requires: txt2mat.m (MatlabCentral File ID: #18430)
%
%   plotFASTOutData() will open a file selector window.  Choose a FAST
%   output file and use the menu to plot desired output variable.
%
%   See also loadFASTOutData

if nargin == 0;
   initialize;
else
   initialize(arg);
end

%------------------

function initialize(arg)
if nargin == 0
    % open file selector and load data
    [A,ffn,nh,SR,hl,fpos] = txt2mat('*.out');
    % return if file selection cancelled
    if isempty(ffn); return; end;
    % parse header lines, assuming labels and units are at end of header
    head = regexp(hl,'[^\n]+','match');
    h.OutList = strtrim(regexp(head{end-1},'[^\t]+','match'));
    h.OutUnits = strtrim(regexp(head{end},'[^\t]+','match'));
    h.OutData = A;
    % extract the filename
    fn =  strtrim(regexp(ffn,'[^\\]+','match'));
    h.fn = fn{end};
    % save data in figure userdata
    set(gcf,'name',h.fn, ...
            'userdata',h);
    
    arg = 2;
end

% get data and select OutData channel to plot
h = get(gcf,'userdata');
mindex = arg;

% plot OutData
clf
id=2;
plot(h.OutData(:,id),h.OutData(:,mindex));
xlabel([h.OutList{id},' ',h.OutUnits{id}]); 
ylabel([h.OutList{mindex},' ',h.OutUnits{mindex}]);

% ui controls for close button and menu
uicontrol(...
   'style','pushbutton', ...
   'units','normalized', ...
   'position',[.70 .94 .12 .06], ...
   'string','close', ...
   'callback','close(gcf)')
uicontrol(...
   'style','popup', ...
   'units','normalized', ...
   'position',[.18 .92 .48 .08], ...
   'string',h.OutList, ...
   'tag','mats', ...
   'fontname','courier', ...
   'fontweight','bold', ...
   'fontsize',10, ...
   'value',mindex, ...
   'callback','plotFASTOutData(get(gco,''value''))');

%------------------
