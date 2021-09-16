function BPESegments(varargin)
%BPESegments  User-interactive determination of BPE segment edge indices
% **********************************************************************
% *           Part of the SNL Wind Turbine Analysis Toolbox            *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% **********************************************************************
%
% 
%

%===== CREDITS & CHANGELOG ================================================
% 2012.05.07  jcb: initial release

%===== INITIALIZATION =====================================================

if ischar(varargin{1})
    % assume filename argument
    nmdfn = varargin{1};
    [a.station, shearweb, active, ansys, a.BladeRotation, blade] = readNuMADinput(nmdfn);
    mp = get(0,'MonitorPositions');
    xpos = (mp(1,3)-800)/2;
    ypos = (mp(1,4)-300)/2;
    app.xy = [xpos ypos 800 300];
    app.job_path = fileparts(nmdfn);
else
    cbo = varargin{1};
    callapp = guidata(cbo);
    app.caller = callapp.fh;  % store figure handle provided by calling script
    app.numadpath = callapp.numadpath;
    app.userpath = callapp.userpath;
    app.job_path = callapp.settings.job_path;
    app.xy = callapp.settings.xy_bpesegments;
    a.station = callapp.station;
    a.BladeRotation = callapp.BladeRotation;
end

N=length(a.station);
BlLength=a.station(end).LocationZ;
app.N = N;

% get BPE segment edges
bpe_sta_ids = fullfile(app.job_path,'bpe_station_ids.txt');
if ~exist(bpe_sta_ids,'file')
    app.bpesta = [];
    app.nbpesta = 0;
else
    fid2 = fopen(bpe_sta_ids,'rt'); % opens basic input file
    if fid2==-1
        errordlg('Could not open "bpe_station_ids.txt"','Error');
        return;
    else
    [app.bpesta app.nbpesta] = fscanf(fid2,'%f,'); % puts values into vector
    app.bpesta = transpose(app.bpesta(:));
    fclose(fid2);
    end
end

chord=[a.station(:).Chord];
z=[a.station(:).LocationZ];
xoff=[a.station(:).Xoffset];

le=-chord.*xoff;
te=-chord.*(xoff-1);
spacing=5; %percent of span
n=1/(spacing*0.01)+1;
guides=linspace(0,BlLength,n);

%===== GUI CONSTRUCTION ===================================================
gui.tag = 'guiBPESegments';
gui.bgcolor = get(0,'defaultUicontrolBackgroundColor');
gui.fig.h = findobj('tag',gui.tag);      % look for tag in current figures
if isempty(gui.fig.h)
    % tag not found, create new
    gui.fig.h = figure('tag',gui.tag,'integerhandle','off','Visible','off');
else
    % tag found, clear figure
    clf(gui.fig.h);
end
set(gui.fig.h,...
    'Name','BPE Segments',...
    'Visible','on',...
    'NextPlot','new',... %causes plot commands to ignore this figure
    'Resize','on',...
    'Units','pixels',...    % pixels | characters
    'Menubar','none',...
    'NumberTitle','off',...
    'Color',gui.bgcolor,...
    'Position',app.xy);
set(gui.fig.h,'CloseRequestFcn',@cb_saveAndClose);
    
rowHeight = 18;
uicontrol(gui.fig.h,...
    'Style','text',...
    'String','Click on stations to select / deselect, or enter station numbers:',...
    'HorizontalAlignment','left',...
    'FontWeight','bold',...
    'Position',[40 80 400 rowHeight]);
str = sprintf('%d, ',app.bpesta);
str(end-1:end) = [];
app.h_textbox = uicontrol(gui.fig.h,...
    'Style','edit',...
    'String',str,...
    'Tag','bpestations',...
    'BackgroundColor','white',...
    'Position',[400 80 350 rowHeight],...
    'Callback',@cb_textbox);
uicontrol(gui.fig.h,...
    'Style','text',...
    'String','Blue triangles are positioned at recommended 5% span spacing.',...
    'HorizontalAlignment','left',...
    'Position',[40 60 400 rowHeight]);
uicontrol(gui.fig.h,...
    'Style','text',...
    'String','Solid green lines are stations selected for BPE analysis.',...
    'HorizontalAlignment','left',...
    'Position',[40 40 400 rowHeight]);
rowHeight = 30;
uicontrol(gui.fig.h,...
    'Style','pushbutton',...
    'String','Save and Close',...
    'Position',[400 40 200 rowHeight],...
    'Callback',@cb_saveAndClose);
uicontrol(gui.fig.h,...
    'Style','pushbutton',...
    'String','Cancel',...
    'Position',[610 40 140 rowHeight],...
    'Callback',@cb_close);    

gui.ax.h = axes('Units','pixels',...
    'Position',[10 100 app.xy(3)-20 app.xy(4)-100]);
set(gui.ax.h,'Units','normalized');
axis equal off;
line([-0.01 1.01]*z(end),[0 0],'Color',gui.bgcolor); % line to create extra space at root & tip
line(guides,ones(n)*1.4*min(le),'Color',[0 0 1],'LineStyle','none','Marker','^');
line(z,le,'Color',[0 0 0]);
line(z,te,'Color',[0 0 0]);
for i=1:N
    app.sta_line(i) = line([z(i) z(i)],[le(i) te(i)]);
    if any(app.bpesta==i)
        set(app.sta_line(i),'Tag','station-bpe','Color',[0 0.8 0],'LineStyle','-','LineWidth',2);
    else
        set(app.sta_line(i),'Tag','station','Color',[0 0 0],'LineStyle',':','LineWidth',1);
    end
    text(z(i),1.8*min(le),num2str(i),'HorizontalAlignment','center','FontSize',6)
end
set(app.sta_line,'ButtonDownFcn',@bdf_selectline);

guidata(gui.fig.h,app);  % store the application data as guidata

end %END GUI CONSTRUCTION (main function)


%===== UTILITY & CALLBACK FUNCTIONS =======================================
function bdf_selectline(cbo,~)
    app = guidata(cbo);
    tag = get(cbo,'Tag');
    
    ksta = find(cbo==app.sta_line);
    
    switch tag
        case 'station-bpe'
            set(cbo,'Tag','station','Color',[0 0 0],'LineStyle',':','LineWidth',1);
            app.bpesta(ksta==app.bpesta) = [];
        case 'station'
            set(cbo,'Tag','station-bpe','Color',[0 0.8 0],'LineStyle','-','LineWidth',2);
            app.bpesta = unique([app.bpesta ksta]);
        otherwise
            return
    end
    
    str = sprintf('%d, ',app.bpesta);
    str(end-1:end) = [];
    set(app.h_textbox,'String',str);
    
    guidata(cbo,app);
end

function cb_textbox(cbo,~)
    app = guidata(cbo);
    
    str = get(cbo,'String');     % get the new input
    pattern = '[\d]*'; % regexp for integers
    match = regexp(str,pattern,'match');
    
    values = zeros(1,numel(match));
    for k=1:numel(match)
        values(k) = str2double(match{k});
    end
    
    outofrange = (values<1) | (values>app.N) | isnan(values);
    values(outofrange) = [];
    app.bpesta = unique(values);
    str = sprintf('%d, ',app.bpesta);
    str(end-1:end) = [];
    set(cbo,'string',str);
    
    set(app.sta_line,'Tag','station','Color',[0 0 0],'LineStyle',':','LineWidth',1);
    set(app.sta_line(app.bpesta),'Tag','station-bpe','Color',[0 0.8 0],'LineStyle','-','LineWidth',2);
    
    guidata(cbo,app);
end

function cb_saveAndClose(cbo,~)
    app = guidata(cbo);  % load application data
    bpe_sta_ids = fullfile(app.job_path,'bpe_station_ids.txt');
    fid2 = fopen(bpe_sta_ids,'wt');
    if fid2==-1
        errordlg('Could not open "bpe_station_ids.txt"','Error');
        return;
    else
        str = sprintf('%d,',app.bpesta);
        str(end) = [];
        try
            fprintf(fid2,'%s  %% NuMAD station numbers to use as edges of BPE segments',str);
        catch ME
            fclose(fid2);
            rethrow(ME);
        end
        fclose(fid2);
    end
    
    cb_close(cbo,[]);
end

function cb_close(cbo,~)
    app = guidata(cbo);  % load application data
    try
        if isfield(app,'caller')
            position = get(gcbf,'Position');
            settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
            if ~isequal(settings.xy_bpesegments,position(1:4))
                settings.xy_bpesegments = position(1:4);
                writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
                callapp = guidata(app.caller);
                callapp.settings.xy_bpesegments = position(1:4);
                guidata(app.caller,callapp);
            end
        end
    catch ME
        h=warndlg('Could not complete close request.');
        uiwait(h);
        closereq;
        rethrow(ME);
    end
    closereq;
end
    
 
  