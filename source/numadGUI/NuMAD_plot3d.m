function NuMAD_plot3d(cbo,~)
%NUMAD_PLOT3D  User interface for Plot3D output options
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_plot3d launches the interface which allows the user to configure
%   Plot3D output options.
%

%===== INITIALIZATION =====================================================
if ~exist('cbo','var')
    errordlg('Script must be called from main NuMAD window.','Programming Error');
    return;
else
    callapp = guidata(cbo);
    app.caller = callapp.fh;  % store figure handle provided by calling script
end

if isfield(callapp,'plot3d');
    app.plot3d = callapp.plot3d;
else
    app.plot3d.n_panels = 100;
    app.plot3d.spacing = 'cosine';
    app.plot3d.newspanloc = [];
    app.plot3d.interpmethod = 'linear';
end

app.debugging = 0;  % with debugging==1, the data structure is sent to the
% workspace after the gui is constructed
%===== GUI CONSTRUCTION ===================================================
app.userpath = callapp.userpath;
app.settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
c=cgui([],'init');  c.fig.xy = app.settings.xy_plot3d;
c=cgui(c,'figure','NuMAD - Plot3D options',[410 275],'guiNuMAD-Plot3D');
set(c.fig.h,'CloseRequestFcn',@cb_close);
app.fh = c.fig.h;

c=cgui(c,'panel','Plot3D options',[10 10 390 260]);
c=cgui(c,'rowHeight','',20);
c=cgui(c,'colEdges','',[10 250 380]);
c=cgui(c,'input','Number of panels per blade surface:',app.plot3d.n_panels,'n_panels',@cb_n_panels);
str = sprintf('%g, ',app.plot3d.breakpoints*100);
str(end-1:end) = [];
c=cgui(c,'input','Breakpoints (% chord):',str,'breakpoints',@cb_breakpoints);
c=cgui(c,'select','Airfoil resampling algorithm:',{'cosine','half-cosine','constant'},'spacing',@cb_spacing);
set(c.ctrl.h,'UserData',{'cosine','half-cosine','constant'});
c=cgui(c,'select','Additional station interpolation method:',{'linear','spline'},'interpmethod',@cb_interpmethod);
set(c.ctrl.h,'UserData',{'linear','spline'});

strspanloc = cell(1,numel(callapp.station));
app.spanloc = zeros(1,numel(callapp.station));
for k=1:numel(callapp.station)
    strspanloc{k} = sprintf('%g',callapp.station(k).LocationZ);
    app.spanloc(k) = callapp.station(k).LocationZ;
end
c=cgui(c,'colEdges','',[10 190]);
c=cgui(c,'listbox','Current station locations (m):',strspanloc,'spanloc','',100);

c=cgui(c,'vpos','previous');
c=cgui(c,'colEdges','',[200 380]);
app.newspanloc = sprintf('%g\n',app.plot3d.newspanloc);
c=cgui(c,'textbox2','Additional station locations (m):',app.newspanloc,'newspanloc',@cb_newspanloc,100);


c=cgui(c,'vpos','value',30);
c=cgui(c,'colEdges','',[10 250]);
c=cgui(c,'push','Save Changes','','',@cb_save);
c=cgui(c,'colEdges','',[250 380]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Discard Changes','','',@cb_close);
% jcb: rather that just using 'close', we should maybe check if anything 
% has changed and confirm with user that they wish to discard changes

% initialize drop-down menus
uictrl = get_uictrl(app.fh,'spacing');
k=find(strcmp(app.plot3d.spacing,uictrl.userdata)==1);
set(uictrl.h,'Value',k);
uictrl = get_uictrl(app.fh,'interpmethod');
k=find(strcmp(app.plot3d.interpmethod,uictrl.userdata)==1);
set(uictrl.h,'Value',k);

guidata(app.fh,app);  % store the application data as guidata

if app.debugging
    assignin('base','app_plot3d',app); 
end

end %END GUI CONSTRUCTION (main function)

%===== UTILITY & CALLBACK FUNCTIONS =======================================
function gui_disable(fh,tags)
    for k=1:numel(tags)
        h = findobj(fh,'tag',tags{k});
        set(h,'Enable','off');
    end
end

function gui_enable(fh,tags)
    for k=1:numel(tags)
        h = findobj(fh,'tag',tags{k});
        set(h,'Enable','on');
    end
end

function gui_invisible(fh,tags)
    for k=1:numel(tags)
        h = findobj(fh,'tag',tags{k});
        set(h,'Visible','off');
    end
end

function gui_visible(fh,tags)
    for k=1:numel(tags)
        h = findobj(fh,'tag',tags{k});
        set(h,'Visible','on');
    end
end

function uictrl = get_uictrl(fh,tag)
    uictrl.h = findobj(fh,'Tag',tag);
    uictrl.string = get(uictrl.h,'String');
    uictrl.value = get(uictrl.h,'Value');
    uictrl.userdata = get(uictrl.h,'UserData');
    if any(strcmp(get(uictrl.h,'Style'),{'popupmenu','listbox'}))
        uictrl.select = uictrl.string(uictrl.value);
    end
end

function cb_n_panels(cbo,~)
    app = guidata(cbo);  % load application data
    val = app.plot3d.n_panels;
    str = get(cbo,'string');     % get the new input
    num = str2double(str);  % try to convert to number
    if isnan(num) || (num<=0)
        % input was not valid => revert to stored value
        str = sprintf('%g',val);
        set(cbo,'string',str);
    else
        app.plot3d.n_panels = floor(num);
        str = sprintf('%g',app.plot3d.n_panels);
        set(cbo,'string',str);
    end
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_plot3d',app); 
    end
end

function cb_breakpoints(cbo,~)
    app = guidata(cbo);  % load application data
    %val = app.plot3d.breakpoints;
    str = get(cbo,'string');     % get the new input
    pattern = '[-+]{0,1}[\d]*[.]*[\d]*[eE]{0,1}[-+]{0,1}[\d]*'; % regexp for a number
    match = regexp(str,pattern,'match');
    
    values = zeros(1,numel(match));
    for k=1:numel(match)
        values(k) = str2double(match{k});
    end
    
    outofrange = (values<=-100) | (values>=100) | isnan(values);
    values(outofrange) = [];
    str = sprintf('%g, ',values);
    str(end-1:end) = [];
    set(cbo,'string',str);
    
    app.plot3d.breakpoints = sort(values/100);
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_plot3d',app); 
    end
end

function cb_spacing(cbo,~)
    app = guidata(cbo);  % load application data
    uictrl.value = get(cbo,'Value');
    uictrl.string = get(cbo,'UserData');  % note: UserData, not String
    uictrl.select = uictrl.string{uictrl.value};
    app.plot3d.spacing = uictrl.select;
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_plot3d',app); 
    end
end

function cb_interpmethod(cbo,~)
    app = guidata(cbo);  % load application data
    uictrl.value = get(cbo,'Value');
    uictrl.string = get(cbo,'UserData');  % note: UserData, not String
    uictrl.select = uictrl.string{uictrl.value};
    app.plot3d.interpmethod = uictrl.select;
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_plot3d',app); 
    end
end

function cb_newspanloc(cbo,~)
    app = guidata(cbo);  % load application data
    str = get(cbo,'String');
    str = [str repmat(' ',size(str,1),1)];  % padd the right side with spaces
    %assignin('base','str',str);
    values = sscanf(transpose(str),'%g');
    outofrange = (values<app.spanloc(1)) | (values>app.spanloc(end));
    for k=1:numel(values)
        if ~isempty(find(values(k)==app.spanloc,1))
            outofrange(k) = 1;
        end
    end
    values(outofrange) = [];
    app.plot3d.newspanloc = sort(values);
    str = sprintf('%g\n',sort(values));
    set(cbo,'String',str(1:end-1)); % update the display
        
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_plot3d',app); 
    end
end

function cb_save(cbo,~)
    app = guidata(cbo);  % load application data
    
    callapp = guidata(app.caller);
    callapp.plot3d = app.plot3d;
    guidata(callapp.fh,callapp);  % store the updated application data
    if callapp.debugging
        assignin('base','app',callapp);
    end
    
    cb_close(cbo,[]);
end

function cb_close(cbo,~)
    try
        app = guidata(cbo);  % load application data
        position = get(app.fh,'Position');
        settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
        if ~isequal(settings.xy_plot3d,position(1:2))
            settings.xy_plot3d = position(1:2);
            writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
        end
    catch ME
        h=warndlg('Could not complete close request.');
        uiwait(h);
        closereq;
        rethrow(ME);
    end
    closereq;
end