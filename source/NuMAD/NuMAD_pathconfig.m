function NuMAD_pathconfig(cbo,~,update_paths)
%NUMAD_PATHCONFIG  User interface for configuring program paths 
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_pathconfig launches the interface which allows the user to set
%   the path to ANSYS, BModes, and PreComp.
%

%===== INITIALIZATION =====================================================
if ~exist('cbo','var')
    errordlg('Script must be called from main NuMAD window.','Programming Error');
    return;
else
    callapp = guidata(cbo);
    app.caller = callapp.fh;  % store figure handle provided by calling script
    app.userpath = callapp.userpath;
end

app.debugging = 0;  % with debugging==1, the data structure is sent to the 
                    % workspace after the gui is constructed
%===== GUI CONSTRUCTION ===================================================
app.settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
app.pathconfig.ansys_path = app.settings.ansys_path;
app.pathconfig.ansys_product = app.settings.ansys_product;
app.pathconfig.bmodes_path = app.settings.bmodes_path;
app.pathconfig.precomp_path = app.settings.precomp_path;

% attempt to read ANSYS product variable file
fn_prodvars = fullfile(app.userpath,'ansys-productvars.txt');
if isequal(2,exist(fn_prodvars,'file'))
    % product variable file found
    fid = fopen(fn_prodvars,'rt');
    if (fid == -1)
        error('Could not open file "%s"',fn_prodvars);
    end
    try
        krow = 1;
        while true
            tline = fgetl(fid);
            if ~ischar(tline), break, end  % reached end-of-file
            if ~isempty(tline) && isequal('%',tline(1)) % found commented line
                continue
            end
            % use regular expression to find Name and Description
            pat = '^\s*(?<prod>\S*)\s+(?<descrip>[^%]*)';
            n = regexp(tline, pat, 'names');
            if isempty(n)
                continue; % continue to next line if regexp fails
            end
            app.ansys_products(krow,:) = {n.prod,strtrim(n.descrip)};
            krow = krow + 1;
        end
        fclose(fid);
    catch ME
        fclose(fid);
        rethrow(ME);
    end
    if isequal(1,krow)
        error('ANSYS Product table is empty. Please fix or delete "%s".',fn_prodvars);
    end
else
    % product variable file not found
    error('Could not find file "%s".\nPlease reopen NuMAD to regenerate.',fn_prodvars);
    % PROGRAMMER - override using following if necessary:
    % app.ansys_products = {'ANE3FL','Multiphysics';
    %                       'ANSYS','Mechanical U'};
end

c=cgui([],'init');  c.fig.xy = app.settings.xy_pathconfig;
c=cgui(c,'figure','NuMAD - Program Paths',[630 300],'guiNuMAD-pathconfig');
set(c.fig.h,'CloseRequestFcn',@cb_close);
c=cgui(c,'panel','Program Paths',[10 10 610 280]);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 200]);
c=cgui(c,'push','Auto-locate ANSYS Path','','',@cb_findansys);
c=cgui(c,'colEdges','',[10 80     500]);
c=cgui(c,'display','ANSYS Path:',app.pathconfig.ansys_path,'ansys_path','');
c=cgui(c,'colEdges','',[10 120     300]);
c=cgui(c,'select','ANSYS Product:',app.ansys_products(:,2),'ansys_product','');
ansysproduct_index = strcmp(app.ansys_products(:,1),app.pathconfig.ansys_product);
if ~any(ansysproduct_index)
    errordlg(sprintf('ANSYS Product Variable "%s" not found. Please add to "%s".',app.settings.ansys_product,fn_prodvars));
    ansysproduct_index(1) = 1;  % just select the first one
end
set(c.ctrl.h,'Value',find(ansysproduct_index==1));
c=cgui(c,'colEdges','',[500 560]);
c=cgui(c,'vpos','previous');
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Browse','','',@cb_browseansys);
c=cgui(c,'colEdges','',[560 600]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Clear','','clear_ansys',@cb_clear);

c=cgui(c,'vpos','value',30);
c=cgui(c,'colEdges','',[10 400]);
c=cgui(c,'push','Save Settings','','',{@cb_save,update_paths});  % note: no @ (already a function handle)
c=cgui(c,'colEdges','',[400 600]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Cancel','','',@cb_close);
% jcb: rather that just using 'close', we should maybe check if anything 
% has changed and confirm with user that they wish to discard changes

c=cgui(c,'panel','Advanced',[20 80 590 100]);
c=cgui(c,'colEdges','',[10 100     470]);
c=cgui(c,'display','BModes Path:',app.pathconfig.bmodes_path,'bmodes_path','');
c=cgui(c,'colEdges','',[480 540]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Browse','','',@cb_browsebmodes);
c=cgui(c,'colEdges','',[540 580]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Clear','','clear_bmodes',@cb_clear);

c=cgui(c,'vpos','next');
c=cgui(c,'colEdges','',[10 100     470]);
c=cgui(c,'display','PreComp Path:',app.pathconfig.precomp_path,'precomp_path','');
c=cgui(c,'colEdges','',[480 540]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Browse','','',@cb_browseprecomp);
c=cgui(c,'colEdges','',[540 580]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Clear','','clear_precomp',@cb_clear);



app.fh = c.fig.h;
guidata(app.fh,app);  % store the application data as guidata

if app.debugging
    assignin('base','app_pathconfig',app); 
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

function cb_findansys(cbo,~)
    % Look for standard ANSYS install location in
    % 'Program Files' and 'Program Files (x86)'.
    % example path: 'C:\Program Files\ANSYS Inc\v130\ansys\bin\winx64\ANSYS130.exe'
    if ~ispc
        warndlg('ANSYS auto-locate available only on windows machines. Please browse for ANSYS manually.','Notification');
        return;
    end
    app = guidata(cbo);  % load application data
    ansys_path = 'C:\Program Files\ANSYS Inc\';
    ff = dir(fullfile(ansys_path,'v*'));
    if isempty(ff)
        ansys_path = 'C:\Program Files (x86)\ANSYS Inc\';
        ff = dir(fullfile(ansys_path,'v*'));
    end
    if ~isempty(ff)  
        versions = {};
        for k=1:numel(ff)
            if ff(k).isdir
                versions(end+1) = {ff(k).name};  %#ok
            end
        end
        sortversions = sort(lower(versions));  % ensure list is sorted
        highestversion = sortversions{end};
        ansys_path = fullfile(ansys_path,highestversion,'ansys','bin');
        platform = {'intel','amd','winx64'};
        for k=1:numel(platform)
            if isdir(fullfile(ansys_path,platform{k}))
                ansys_path = fullfile(ansys_path,platform{k});
                break;
            end
        end
        executable = sprintf('ANSYS%s.exe',highestversion(2:end));
        ansys_path = fullfile(ansys_path,executable);
        
        % attempt to read product variable
        ansys_product = '';
        launcher_profile = fullfile(getenv('APPDATA'),'Ansys',highestversion,'launcher','profiles.xml');
        if exist(launcher_profile,'file')
            fid = fopen(launcher_profile,'rt');
            if (fid == -1)
                fprintf('Could not open file "%s"\n',launcher_profile);
            else
                filecontents = fread(fid,inf,'uint8=>char')';
                fclose(fid);
                % try to extract the last-used ansys product
                n=regexp(filecontents,'lic="(?<lic>[^"]*)"','names');
                if ~isempty(n.lic)
                    ansys_product = n.lic;
                end
            end
        end
        
        if exist(ansys_path,'file')
            app.pathconfig.ansys_path = ansys_path;
            if ~isempty(ansys_product)
                ansysproduct_index = strcmpi(app.ansys_products(:,1),ansys_product);
                if any(ansysproduct_index)
                    app.pathconfig.ansys_product = ansys_product;
                    uictrl = get_uictrl(gcbf,'ansys_product');
                    set(uictrl.h,'Value',find(ansysproduct_index==1));
                end
            end
            guidata(app.fh,app);  % store the updated guidata
            if app.debugging
                assignin('base','app_pathconfig',app);
            end
            uictrl = get_uictrl(gcbf,'ansys_path');
            set(uictrl.h,'String',ansys_path);
            helpdlg(sprintf('ANSYS found!\n%s',ansys_path),'Notification');
            return;  % this RETURN skips past the helpdlg below
        end
    end
    warndlg('Could not auto-locate ANSYS. Please browse for ANSYS manually.','Notification');    
end

function cb_browseansys(cbo,~)
%     helpstr = sprintf('Please select the ANSYS executable.\n%s', ...
%     'example: C:\Program Files\ANSYS Inc\v130\ansys\bin\winx64\ANSYS130.exe');
%     h = helpdlg(helpstr,'Notification');
%     waitfor(h);
    if isdir('C:\Program Files\ANSYS Inc\')
        startpath = 'C:\Program Files\ANSYS Inc\';
    elseif isdir('C:\Program Files\')
        startpath = 'C:\Program Files\';
    elseif isdir('C:\');
        startpath = 'C:\';
    else
        startpath = '';
    end
    [fn pn] = uigetfile( ...
        {'ANSYS*.exe', 'ANSYS executable'; '*.*', 'All Files (*.*)'}, ...
        'Example: C:\Program Files\ANSYS Inc\v130\ansys\bin\winx64\ANSYS130.exe',startpath);
    if isequal(fn,0) || isequal(pn,0)
        return
    else
        ansys_path = fullfile(pn,fn);
        app = guidata(cbo);  % load application data
        app.pathconfig.ansys_path = ansys_path;
        guidata(app.fh,app);  % store the updated guidata
        if app.debugging
            assignin('base','app_pathconfig',app);
        end
        uictrl = get_uictrl(gcbf,'ansys_path');
        set(uictrl.h,'String',ansys_path);
    end
end

function cb_browsebmodes(cbo,~)
    if isdir('C:\');
        startpath = 'C:\';
    else
        startpath = '';
    end
    [fn pn] = uigetfile( ...
        {'*.exe', 'Programs (*.exe)'; '*.*', 'All Files (*.*)'}, ...
        'Select BModes executable',startpath);
    if isequal(fn,0) || isequal(pn,0)
        return
    else
        bmodes_path = fullfile(pn,fn);
        app = guidata(cbo);  % load application data
        app.pathconfig.bmodes_path = bmodes_path;
        guidata(app.fh,app);  % store the updated guidata
        if app.debugging
            assignin('base','app_pathconfig',app);
        end
        uictrl = get_uictrl(gcbf,'bmodes_path');
        set(uictrl.h,'String',bmodes_path);
    end
end

function cb_browseprecomp(cbo,~)
    if isdir('C:\');
        startpath = 'C:\';
    else
        startpath = '';
    end
    [fn pn] = uigetfile( ...
        {'*.exe', 'Programs (*.exe)'; '*.*', 'All Files (*.*)'}, ...
        'Select PreComp executable',startpath);
    if isequal(fn,0) || isequal(pn,0)
        return
    else
        precomp_path = fullfile(pn,fn);
        app = guidata(cbo);  % load application data
        app.pathconfig.precomp_path = precomp_path;
        guidata(app.fh,app);  % store the updated guidata
        if app.debugging
            assignin('base','app_pathconfig',app);
        end
        uictrl = get_uictrl(gcbf,'precomp_path');
        set(uictrl.h,'String',precomp_path);
    end
end

function cb_clear(cbo,~)
    app = guidata(cbo);  % load application data
    tag = get(cbo,'tag');
    switch tag
        case 'clear_ansys'
            app.pathconfig.ansys_path = '';
            uictrl = get_uictrl(gcbf,'ansys_path');
            set(uictrl.h,'String','');
        case 'clear_bmodes'
            app.pathconfig.bmodes_path = '';
            uictrl = get_uictrl(gcbf,'bmodes_path');
            set(uictrl.h,'String','');
        case 'clear_precomp'
            app.pathconfig.precomp_path = '';
            uictrl = get_uictrl(gcbf,'precomp_path');
            set(uictrl.h,'String','');
    end
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_pathconfig',app);
    end
end

function cb_save(cbo,~,update_paths)
    app = guidata(cbo);  % load application data
    uictrl = get_uictrl(gcbf,'ansys_product');
    app.pathconfig.ansys_product = app.ansys_products{uictrl.value,1};
    update_paths(app.caller,app.pathconfig);  % send data to caller
    cb_close(cbo,[]);
end

function cb_close(cbo,~)
    app = guidata(cbo);  % load application data
    position = get(app.fh,'Position');
    settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
    if ~isequal(settings.xy_pathconfig, position(1:2))
        settings.xy_pathconfig = position(1:2);
        writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
    end
    closereq;
end