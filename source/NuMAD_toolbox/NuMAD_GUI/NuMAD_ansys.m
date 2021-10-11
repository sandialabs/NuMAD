function NuMAD_ansys(cbo,~,update_ansys)
%NUMAD_ANSYS  User interface to ANSYS output options
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_ansys launches the ANSYS output options user interface.
%   Usage: uimenu('label','Output Options','callback',{@NuMAD_ansys,@update_ansys});
%   where update_ansys is a function in NuMAD_main

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
app.ansys = callapp.ansys;
app.ansys.path = app.settings.ansys_path;
app.ansys.product = app.settings.ansys_product;

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

c=cgui([],'init');  c.fig.xy = app.settings.xy_ansys;
c=cgui(c,'figure','NuMAD - ANSYS options',[410 500],'guiNuMAD-ANSYS');
set(c.fig.h,'CloseRequestFcn',@cb_close);
c=cgui(c,'panel','ANSYS options',[10 10 390 480]);
c=cgui(c,'rowHeight','',30);
c=cgui(c,'colEdges','',[10 80     380]);
c=cgui(c,'display','ANSYS Path:',app.settings.ansys_path,'ansys_path','');
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 200]);
c=cgui(c,'push','Auto-Locate ANSYS','','',@cb_findansys);
c=cgui(c,'colEdges','',[200 380]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Browse for ANSYS','','',@cb_browseansys);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 140     380]);
c=cgui(c,'select','ANSYS Product:',app.ansys_products(:,2),'ansys_product','');
ansysproduct_index = strcmp(app.ansys_products(:,1),app.settings.ansys_product);
if ~any(ansysproduct_index)
    errordlg(sprintf('ANSYS Product Variable "%s" not found. Please add to "%s".',app.settings.ansys_product,fn_prodvars));
    ansysproduct_index(1) = 1;  % just select the first one
end
set(c.ctrl.h,'Value',find(ansysproduct_index==1));
c=cgui(c,'vpos','next');
c=cgui(c,'select','Boundary Conditions:',{'Cantilevered','Free-Free'},'BoundaryCondition',@cb_select);
set(c.ctrl.h,'UserData',{'cantilevered','freefree'});
c=cgui(c,'select','Element System:',{'8-Node Structural Shell (SHELL281)','4-Node Structural Shell (SHELL181)',...
    'Offset Layered Shells (SHELL99)','Nonlinear Offset Layered Shells (SHELL91)'},'ElementSystem',@cb_select);
set(c.ctrl.h,'UserData',{'281','181','99','91'});
c=cgui(c,'select','Mult. Layer Behavior:',{'Produce distinct layers','Multiply layer thickness'},'MultipleLayerBehavior',@cb_select);
set(c.ctrl.h,'UserData',{'distinct','multiply'});
c=cgui(c,'select','Mesh Density Method:',{'Smart Mesh','Element Size'},'meshing',@cb_meshing);
set(c.ctrl.h,'UserData',{'smartsize','elementsize'});
c=cgui(c,'colEdges','',[   140 290 380]);
smartsize = cellstr(num2str(transpose(1:10)));
smartsize([1 10]) = [{' 1 (fine)'};{'10 (coarse)'}];
c=cgui(c,'select','Smart Mesh Parameter:',smartsize,'smartsize',@cb_smartsize);
set(c.ctrl.h,'Value',app.ansys.smartsize);
set(findobj(c.fig.h,'string','Smart Mesh Parameter:'),'tag','smartsize-label');  %jcb: recode this later with new cgui approach
c=cgui(c,'vpos','previous');
c=cgui(c,'input','Element Edge Length:',app.ansys.elementsize,'elementsize',@cb_elementsize);
set(findobj(c.fig.h,'string','Element Edge Length:'),'tag','elementsize-label');  %jcb: recode this later with new cgui approach

tableColumns = {'Failure Criteria Post-Processing', 'Activate'};
tableData = {'Maximum strain criterion (EMAX)','no';...
             'Maximum stress criterion (SMAX)','no';...
             'Tsai-Wu strength index (TWSI)','no';...
             'Inverse of Tsai-Wu strength ratio index (TWSR)','no';...
             'Hashin fiber failure criterion (HFIB)','no';...
             'Hashin matrix failure criterion (HMAT)','no';...
             'Puck fiber failure criterion (PFIB)','no';...
             'Puck inter-fiber (matrix) failure criterion (PMAT)','no';...
             'LaRc03 fiber failure criterion (L3FB)','no';...
             'LaRc04 matrix failure criterion (L3MT)','no';...
             'LaRc04 fiber failure criterion (L4FB)','no';...
             'LaRc04 matrix failure criterion (L4MT)','no';...
             'User-defined failure criterion (USR1)','no';...
             'User-defined failure criterion (USR2)','no';...
             'User-defined failure criterion (USR3)','no';...
             'User-defined failure criterion (USR4)','no';...
             'User-defined failure criterion (USR5)','no';...
             'User-defined failure criterion (USR6)','no';...
             'User-defined failure criterion (USR7)','no';...
             'User-defined failure criterion (USR8)','no';...
             'User-defined failure criterion (USR9)','no'};
if length(tableData)~=length(app.ansys.FailureCriteria)
    error('failure criteria table length mismatch');
end
for kfc = 1:length(app.ansys.FailureCriteria)
    if app.ansys.FailureCriteria{kfc,2}
        tableData{kfc,2} = 'yes';
    end
end
         
app.th=uitable('Parent',c.panel.h,'Position',[10 120 365 100],...
    'ColumnWidth',{255 55},'ColumnEditable',[false true],...
    'Data',tableData,'ColumnName',tableColumns);
% restrict set of choices for material column
set(app.th,'ColumnFormat',{[] {'no','yes'}}); 
% install a callback for data verification
%set(app.th,'CellEditCallback',@cb_layertable);

c=cgui(c,'vpos','next');
c=cgui(c,'vpos','next');
c=cgui(c,'vpos','next');
c=cgui(c,'vpos','next');
c=cgui(c,'vpos','next');
c=cgui(c,'colEdges','',[10 140     380]);
c=cgui(c,'label','Output options:');
c=cgui(c,'vpos','previous');
c=cgui(c,'colEdges','',[   140     280]);
c=cgui(c,'checkbox','Input file (shell7.src)',app.ansys.shell7gen,'shell7gen',@cb_checkbox);
c=cgui(c,'checkbox','ANSYS database (.db)',app.ansys.dbgen,'dbgen',@cb_checkbox);

c=cgui(c,'vpos','next');
c=cgui(c,'colEdges','',[10 250]);
c=cgui(c,'push','Save Changes','','',{@cb_save,update_ansys});  % note: no @ (already a function handle)
c=cgui(c,'colEdges','',[250 380]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Discard Changes','','',@cb_close);
% jcb: rather that just using 'close', we should maybe check if anything 
% has changed and confirm with user that they wish to discard changes

app.fh = c.fig.h;
guidata(app.fh,app);  % store the application data as guidata

% initialize drop-down menus
uictrl = get_uictrl(app.fh,'BoundaryCondition');
k=find(strcmp(app.ansys.BoundaryCondition,uictrl.userdata)==1);
set(uictrl.h,'Value',k);
uictrl = get_uictrl(app.fh,'ElementSystem');
k=find(strcmp(app.ansys.ElementSystem,uictrl.userdata)==1);
set(uictrl.h,'Value',k);
uictrl = get_uictrl(app.fh,'MultipleLayerBehavior');
k=find(strcmp(app.ansys.MultipleLayerBehavior,uictrl.userdata)==1);
set(uictrl.h,'Value',k);
uictrl = get_uictrl(app.fh,'meshing');
k=find(strcmp(app.ansys.meshing,uictrl.userdata)==1);
set(uictrl.h,'Value',k);
cb_meshing(uictrl.h,[]);  % enable/disable other controls based on meshing option

%gui_disable(app.fh,{'shell7gen','dbgen'});

if app.debugging
    assignin('base','app_ansys',app); 
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
            app.ansys.path = ansys_path;
            if ~isempty(ansys_product)
                ansysproduct_index = strcmpi(app.ansys_products(:,1),ansys_product);
                if any(ansysproduct_index)
                    app.ansys.product = ansys_product;
                    uictrl = get_uictrl(gcbf,'ansys_product');
                    set(uictrl.h,'Value',find(ansysproduct_index==1));
                end
            end
            guidata(app.fh,app);  % store the updated guidata
            if app.debugging
                assignin('base','app_ansys',app);
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
        app.ansys.path = ansys_path;
        guidata(app.fh,app);  % store the updated guidata
        if app.debugging
            assignin('base','app_ansys',app);
        end
        uictrl = get_uictrl(gcbf,'ansys_path');
        set(uictrl.h,'String',ansys_path);
    end
end

function cb_select(cbo,~)
    app = guidata(cbo);  % load application data
    tag = get(cbo,'Tag');
    uictrl.value = get(cbo,'Value');
    uictrl.string = get(cbo,'UserData');  % note: UserData, not String
    uictrl.select = uictrl.string{uictrl.value};
    app.ansys.(tag) = uictrl.select;
    
    if isequal('ElementSystem',tag)
        switch uictrl.select
            case {'91','99'}
                % write_shell7 always use the 'multiply' behavior for shell91/99
                h = findobj(app.fh,'Tag','MultipleLayerBehavior');
                set(h,'Value',2,'Enable','off'); % 2 -> 'multiply' option
            otherwise
                % reverting back to shell281/181 uses last stored value
                h = findobj(app.fh,'Tag','MultipleLayerBehavior');
                set(h,'Enable','on');
                switch app.ansys.MultipleLayerBehavior
                    case 'distinct'
                        set(h,'Value',1);  
                    case 'multiply'
                        set(h,'Value',2);
                end
        end
    end
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_ansys',app); 
    end
end

function cb_meshing(cbo,~)
    app = guidata(cbo);  % load application data
    uictrl.value = get(cbo,'Value');
    uictrl.string = get(cbo,'UserData');  % note: UserData, not String
    uictrl.select = uictrl.string{uictrl.value};
    app.ansys.meshing = uictrl.select;
    
    switch uictrl.select
        case 'smartsize'
            gui_visible(app.fh,{'smartsize','smartsize-label'});
            gui_invisible(app.fh,{'elementsize','elementsize-label'});
        case 'elementsize'
            gui_invisible(app.fh,{'smartsize','smartsize-label'});
            gui_visible(app.fh,{'elementsize','elementsize-label'});
    end
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_ansys',app); 
    end
end

function cb_smartsize(cbo,~)
    app = guidata(cbo);  % load application data
    app.ansys.smartsize = get(cbo,'Value');
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_ansys',app); 
    end
end

function cb_elementsize(cbo,~)
    app = guidata(cbo);  % load application data
    val = app.ansys.elementsize;
    str = get(cbo,'string');     % get the new input
    num = str2double(str);  % try to convert to number
    if isnan(num) || (num<=0)
        % input was not valid => revert to stored value
        str = sprintf('%g',val);
        set(cbo,'string',str);
    else
        app.ansys.elementsize = num;
    end
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_ansys',app); 
    end
end

function cb_checkbox(cbo,~)
    app = guidata(cbo);  % load application data
    tag = get(cbo,'Tag');
    value = get(cbo,'Value');

    app.ansys.(tag) = value;
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_ansys',app); 
    end
end

function cb_save(cbo,~,update_ansys)
    app = guidata(cbo);  % load application data
    if ~app.ansys.shell7gen && ~app.ansys.dbgen
        helpdlg('Please select at least one output option.','Notification');
        return;
    end
    % read changes to ansys product
    uictrl = get_uictrl(gcbf,'ansys_product');
    app.ansys.product = app.ansys_products{uictrl.value,1};
    % read changes to failure criteria table
    tableData = get(app.th,'Data');
    for kfc = 1:length(tableData)
        app.ansys.FailureCriteria{kfc,2} = isequal('yes',tableData{kfc,2});
    end
    update_ansys(app.caller,app.ansys);  % send data to caller
    cb_close(cbo,[]);
end

function cb_close(cbo,~)
    app = guidata(cbo);  % load application data
    position = get(app.fh,'Position');
    settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
    if ~isequal(settings.xy_ansys,position(1:2))
        settings.xy_ansys = position(1:2);
        writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
    end
    closereq;
end