function NuMAD_activelist(cbo,~,update_activelist)
%NUMAD_ACTIVELIST  NuMAD material palette user interface
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_activelist launches the interface which allows the user to
%   select a subset of the material database to use while assigning
%   materials to blade sections.
%   Usage: uimenu('label','Material Palette','callback',...
%                 {@NuMAD_activelist,@update_activelist});
%   where update_activelist is a function in NuMAD_main


if ~exist('cbo','var')
    errordlg('Script must be called from main NuMAD window.','Programming Error');
    return;
else
    callapp = guidata(cbo);
    app.caller = callapp.fh;  % store figure handle provided by calling script
    app.numadpath = callapp.numadpath;
    app.userpath = callapp.userpath;
end

app.debugging = 0;  % with debugging==1, the data structure is sent to the 
                    % workspace after the gui is constructed
%===== GUI CONSTRUCTION ===================================================
%units = unitconst();

app.settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
c=cgui([],'init');  c.fig.xy = app.settings.xy_activelist;
c=cgui(c,'figure','NuMAD - select material palette (from MASTER)',[780 260],'guiNuMAD-activelist');
set(c.fig.h,'CloseRequestFcn',@cb_close);


if exist('callapp','var') && ~isempty(callapp.settings.job_name)
    % use the local material database
    app.matdb_path = fullfile(callapp.settings.job_path,'MatDBsi.txt');
    set(c.fig.h,'Name','NuMAD - select material palette (from LOCAL)');
else
    % if no job name, use the master material database
    app.matdb_path = fullfile(app.numadpath,'MatDBsi.txt');
end
app.matdb = readMatDB(app.matdb_path);
for k=1:numel(app.matdb)
    app.matlist{k} = app.matdb(k).name;
    mattype{k} = app.matdb(k).type;
end
app.isotropic =  strcmp('isotropic',mattype);
app.orthotropic = strcmp('orthotropic',mattype);
app.composite = strcmp('composite',mattype);
app.selected = '';
app.active = callapp.active;

c=cgui(c,'panel','',[10 10 760 240]);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'vpos','value',230);
c=cgui(c,'colEdges','',[10 250]);
c=cgui(c,'listbox','All Composite Materials',sort(app.matlist(app.composite)),'Composites','',180);

c=cgui(c,'vpos','value',180);
c=cgui(c,'colEdges','',[250 350]);
c=cgui(c,'push','Add =>','','',@cb_add);
c=cgui(c,'vpos','next');
c=cgui(c,'push','Remove <=','','',@cb_remove);

c=cgui(c,'vpos','value',230);
c=cgui(c,'colEdges','',[350 600]);
c=cgui(c,'listbox','Current Material Palette',app.active.list,'ActiveList','',180);

c=cgui(c,'vpos','value',180);
c=cgui(c,'colEdges','',[600 700]);
c=cgui(c,'push','Save Changes','','',{@cb_save,update_activelist});

% app.colors = {[255   0   0]/255;  % red
%               [  0 102 255]/255;  % blue
%               [  0 204   0]/255;  % green
%               [204 102   0]/255;  % brown
%               [255 255   0]/255;  % yellow
%               [  0 153 255]/255;  % sky blue
%               [255 153   0]/255;  % orange
%               [  0 153   0]/255;  % dark green
%               [204 102 255]/255;  % purple
%               [204 153   0]/255;  % light brown
%               [  0 255 255]/255;  % cyan
%               [255  51   0]/255;  % light red
%               [  0 255   0]/255;  % bright green
%               [204   0   0]/255;  % dark red
%               [255 204 102]/255;  % tan
%               [  0 102   0]/255;  % forest green
%               [255   0 255]/255;  % magenta
%               [102 102 102]/255;  % gray - 40%
%               [255   0 102]/255;  % name? (pinkish red)
%               [102 153   0]/255}; % olive green

app.fh = c.fig.h;
guidata(app.fh,app);  % store the application data as guidata

if app.debugging
    assignin('base','app_act',app); 
end

end %END GUI CONSTRUCTION (main function)


%===== UTILITY & CALLBACK FUNCTIONS =======================================
function uictrl = get_uictrl(fh,tag)
    uictrl.h = findobj(fh,'Tag',tag);
    uictrl.string = get(uictrl.h,'String');
    uictrl.value = get(uictrl.h,'Value');
    uictrl.userdata = get(uictrl.h,'UserData');
    if any(strcmp(get(uictrl.h,'Style'),{'popupmenu','listbox'}))
        uictrl.select = uictrl.string(uictrl.value);
    end
end

function cb_add(cbo,~)
    app = guidata(cbo);  % load application data
    
    % get the menulist and selection of the listbox
    
    comp = get_uictrl(app.fh,'Composites');
    if numel(comp.select)>0
        newlist = [app.active.list; comp.select];
        app.active.list = unique(newlist);
%         kclr = rem((1:numel(app.active.list))-1,numel(app.colors))+1;
%         app.active.color = app.colors(kclr);
    else
        return;
    end
    
    obj = findobj(app.fh,'Tag','ActiveList');
    set(obj,'String',app.active.list);
    
    guidata(app.fh,app);  % store the application data as guidata
    if app.debugging
        assignin('base','app_act',app);
    end
end

function cb_remove(cbo,~)
    app = guidata(cbo);  % load application data
    
    % get the menulist and selection of the listbox
    
    actlist = get_uictrl(app.fh,'ActiveList');
    newlist = actlist.string;
    newlist(actlist.value) = [];
    app.active.list = unique([newlist; '**UNSPECIFIED**']);
%     kclr = rem((1:numel(app.active.list))-1,numel(app.colors))+1;
%     app.active.color = app.colors(kclr);
    set(actlist.h,'String',app.active.list,'Value',[]);
    
    guidata(app.fh,app);  % store the application data as guidata
    if app.debugging
        assignin('base','app_act',app);
    end
end

function cb_save(cbo,~,update_activelist)
    app = guidata(cbo);  % load application data
    update_activelist(app.caller,app.active);  % send data to caller
    cb_close(cbo,[]);
end

function cb_close(cbo,~)
    app = guidata(cbo);  % load application data
    position = get(app.fh,'Position');
    settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
    if ~isequal(settings.xy_activelist,position(1:2))
        settings.xy_activelist = position(1:2);
        writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
    end
    closereq;
end
