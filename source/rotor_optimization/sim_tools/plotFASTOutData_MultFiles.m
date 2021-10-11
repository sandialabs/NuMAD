function plotFASTOutData_MultFiles(arg)
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
    % SPECIFY FILENAMES MANUALLY
    [A,ffn,nh,SR,hl,fpos] = txt2mat('C:\Users\brresor\Documents\SystemModels\SNL_3MW_CRP_90\IECSweep_SNL3MW.out');
    B = txt2mat('C:\Users\brresor\Documents\SystemModels\SNL_3MW_CRP_80\IECSweep_SNL3MW.out');
    C = txt2mat('NoCup_BRR_15MW_TorStudy.out');
    D = txt2mat('Test13_ADAMS.plt');
    if (isempty(A) || isempty(B) || isempty(C) || isempty(D))
        return
    end
    % parse header lines, assuming labels and units are at end of header
    head = regexp(hl,'[^\n]+','match');
    h.OutList = strtrim(regexp(head{end-1},'[^\t]+','match'));
    h.OutUnits = strtrim(regexp(head{end},'[^\t]+','match'));
    h.A = A;
    h.B = B;
    h.C = C;
    h.D = D;
    % save data in figure userdata
    set(gcf,'userdata',h);
    
    arg = 2;
end

% get data and select OutData channel to plot
h = get(gcf,'userdata');
mindex = arg;

% plot OutData
clf
plot(h.A(:,1),h.A(:,mindex),...
     h.B(:,1),h.B(:,mindex),...
     h.C(:,1),h.C(:,mindex),...
     h.D(:,1),h.D(:,mindex));
xlabel([h.OutList{1},' ',h.OutUnits{1}]); 
ylabel([h.OutList{mindex},' ',h.OutUnits{mindex}]);
legend('FAST','CurveFAST TorsOff','CurveFAST TorsOn','ADAMS');

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
   'callback','plotFASTOutData_MultFiles(get(gco,''value''))');

%------------------
