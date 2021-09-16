function NuMAD_materials(cbo,~)
%NUMAD_MATERIALS  User interface for materials database
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_isotropic launches the user interface for the materials database.
%

%===== INITIALIZATION =====================================================
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
c=cgui([],'init');  c.fig.xy = app.settings.xy_materials;
c=cgui(c,'figure','NuMAD - material database (MASTER)',[780 290],'guiNuMAD-materials');
set(c.fig.h,'CloseRequestFcn',@cb_close);

if exist('callapp','var') && ~isempty(callapp.settings.job_name)
    % use the local material database
    app.matdb_path = fullfile(callapp.settings.job_path,'MatDBsi.txt');
    set(c.fig.h,'Name','NuMAD - material database (LOCAL)');
else
    % if no job name, use the master material database
    app.matdb_path = fullfile(app.numadpath,'MatDBsi.txt');
end
app.matdb = readMatDB(app.matdb_path);
% success = copyfile('MatDBsi.txt','MatDBsi.bak');
% if ~success
%     warndlg('Unable to backup material database','Warning');
% end
for k=1:numel(app.matdb)
    app.matlist{k} = app.matdb(k).name;
    mattype{k} = app.matdb(k).type;
end
app.isotropic =  strcmp('isotropic',mattype);
app.orthotropic = strcmp('orthotropic',mattype);
app.composite = strcmp('composite',mattype);
app.selected = '';

f=uimenu('label','New');
  uimenu(f,'label','Isotropic','callback',{@cb_new,'isotropic'});
  uimenu(f,'label','Ortho/Layer','callback',{@cb_new,'orthotropic'});
  uimenu(f,'label','Composite','callback',{@cb_new,'composite'});
  uimenu(f,'label','Duplicate Selected','callback',@cb_duplicate,'separator','on');
f=uimenu('label','Edit');
  uimenu(f,'label','Modify Selected','callback',@cb_modify);
  uimenu(f,'label','Delete Selected','callback',@cb_delete);

c=cgui(c,'panel','',[10 10 760 240]);
c=cgui(c,'vpos','value',230);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 250  400]);
c=cgui(c,'listbox','Isotropic',app.matlist(app.isotropic),'iso',@cb_matlist,180);
c=cgui(c,'colEdges','',[250 500]);
c=cgui(c,'vpos','value',230);
c=cgui(c,'listbox','Orthotropic/Layer',app.matlist(app.orthotropic),'ortho',@cb_matlist,180);
c=cgui(c,'colEdges','',[500 750]);
c=cgui(c,'vpos','value',230);
c=cgui(c,'listbox','Composite',app.matlist(app.composite),'comp',@cb_matlist,180);

app.fh = c.fig.h;
guidata(app.fh,app);  % store the application data as guidata

if app.debugging
    assignin('base','app_mat',app); 
end

end %END GUI CONSTRUCTION (main function)


%===== UTILITY & CALLBACK FUNCTIONS =======================================
function refresh_matlistboxes(fh)
    app = guidata(fh);  % load application data
    
    [app.matlist mattype] = deal(cell(1,numel(app.matdb)));  %initialize
    for k=1:numel(app.matdb)
        app.matlist{k} = app.matdb(k).name;
        mattype{k} = app.matdb(k).type;
    end
    app.isotropic =  strcmp('isotropic',mattype);
    app.orthotropic = strcmp('orthotropic',mattype);
    app.composite = strcmp('composite',mattype);
    app.selected = '';
    
    % find the three listboxes
    iso = findobj(app.fh,'Tag','iso');
    ortho = findobj(app.fh,'Tag','ortho');
    comp = findobj(app.fh,'Tag','comp');
    
    set(iso,'String',app.matlist(app.isotropic));
    set(ortho,'String',app.matlist(app.orthotropic));
    set(comp,'String',app.matlist(app.composite));
    set([iso,ortho,comp],'Value',[]);
    
    guidata(app.fh,app);  % store the application data as guidata
    if app.debugging
        assignin('base','app_mat',app);
    end
end

function cb_matlist(cbo,~)
% this callback responds to user interaction with the three material
% listboxes (isotropic, orthotropic, composite)
    app = guidata(cbo);  % load application data
    % find the three listboxes
    iso = findobj(app.fh,'Tag','iso');
    ortho = findobj(app.fh,'Tag','ortho');
    comp = findobj(app.fh,'Tag','comp');
    % these IF statements create the gui behavior of allowing
    % selection in only one of the three listboxes
    if cbo~=iso, set(iso,'Value',[]); end
    if cbo~=ortho, set(ortho,'Value',[]); end    
    if cbo~=comp, set(comp,'Value',[]); end
    
    % get the menulist and selection of the current listbox
    menulist = get(cbo,'string');
    selection = get(cbo,'value');
    if numel(selection)==1
        % if only one item is selected, remember it
        app.selected = menulist{selection};
    else
        % else, do not remember selection
        app.selected = '';
    end
    
    guidata(app.fh,app);  % store the application data as guidata
    if app.debugging
        assignin('base','app_mat',app);
    end
    
    % these lines allow the user to double-click on a material name
    % to open it for modification
    click = get(gcbf,'SelectionType');
    switch click
        case 'open'
        cb_modify(cbo,[]);    
    end
end

function cb_new(cbo,~,type)
% this callback creates a new material and opens it
% with the appropriate script for modification
    app = guidata(cbo);  % load application data
    
    
    % jcb: this structure must match the one in 'readMatDB.m'
    newmat = deal(struct('type',[],'name',[],'reference',[],...
    'dens',[],'nuxy',[],'ex',[],'ey',[],'ez',[],'gxy',[],'gyz',[],'gxz',[],...
    'prxy',[],'pryz',[],'prxz',[],'xten',[],'xcmp',[],'yten',[],'ycmp',[],...
    'zten',[],'zcmp',[],'xy',[],'yz',[],'xz',[],'xycp',[],'yzcp',[],'xzcp',[],...
    'xzit',[],'xzic',[],'yzit',[],'yzic',[],'g1g2',[],'etal',[],'etat',[],...
    'alp0',[],'thicknessType',[],'uniqueLayers',[],'symmetryType',[],'layer',[]));
    newmat.type = type;
    % each script should look for '__NewMatFlag' and then populate the
    % appropriate fields
    newmat.name = '__NewMatFlag';
    newmat.reference = '';
    
    switch type
        % call the appropriate gui, providing the: 
        %   1) material data
        %   2) the fh of the main material window
        %   3) and the handle of the function which saves 
        %      the data back into the database
        case 'isotropic'
            NuMAD_isotropic(newmat,app.fh,@cb_matsave);
        case 'orthotropic'
            NuMAD_orthotropic(newmat,app.fh,@cb_matsave);
        case 'composite'
            NuMAD_composite(newmat,app.fh,@cb_matsave);
    end
end

function cb_modify(cbo,~)
% this callback opens a material with the appropriate script for
% modification
    app = guidata(cbo);  % load application data
    
    if isempty(app.selected)
        helpdlg('Select a single material to modify.','Operation Not Permitted');
        return;
    end
    
    % find the selection in the material list
    k = strcmp(app.selected,app.matlist);
    % load the material data from the database
    mat = app.matdb(k);
        
    switch mat.type
        % call the appropriate gui, providing the: 
        %   1) material data
        %   2) the fh of the main material window
        %   3) and the handle of the function which saves 
        %      the data back into the database
        case 'isotropic'
            NuMAD_isotropic(mat,app.fh,@cb_matsave);
        case 'orthotropic'
            NuMAD_orthotropic(mat,app.fh,@cb_matsave);
        case 'composite'
            NuMAD_composite(mat,app.fh,@cb_matsave);
    end
end

function cb_duplicate(cbo,~)
% this callback saves a material definition with a new name
    app = guidata(cbo);  % load application data
    
    if isempty(app.selected)
        helpdlg('Select a single material to duplicate.','Operation Not Permitted');
        return;
    end
    
    % find the selection in the material list
    k = strcmp(app.selected,app.matlist);
    % load the material data from the database
    mat = app.matdb(k);

    defaultanswer = {mat.name};
    while true
        % note: inputdlg returns a cell string
        answer = inputdlg('Duplicate the selected material and save with the following name:','Duplicate',1,defaultanswer);
        
        if isempty(answer)
            return;  % the user canceled the operation
        end
        
        % check if the new name is already in the material list
        k2 = strcmp(answer,app.matlist);
        if any(k2)
            h=warndlg('That material name already exists.','Operation Not Permitted');
            waitfor(h);  % wait until the user closes the warning dialog
            defaultanswer = answer;  % insert the last answer into the inputdlg
        else
            break  % the new name is acceptable
        end
    end

    mat.name = answer{1};  % note: inputdlg returns a cell string
    
    % append the new material to end of matdb
    k = numel(app.matdb) + 1;  % note k is now numerical index rather than logical
    app.matdb(k) = mat;
    
    % write out the changes to the material database file
    writeMatDB(app.matdb,app.matdb_path);
    refresh_matdb(app)
    
    guidata(app.fh,app);  % store the application data as guidata
    if app.debugging
        assignin('base','app_mat',app);
    end
    
    refresh_matlistboxes(app.fh);
end

function cb_delete(cbo,~)
% this function deletes the currently selected material
    app = guidata(cbo);  % load application data
    
    if isempty(app.selected)
        helpdlg('Select a single material to delete.','Operation Not Permitted');
        return;
    end
    
    % find the selection in the material list
    k = strcmp(app.selected,app.matlist);
    % load the material data from the database
    matname = app.matdb(k).name;

    % determine which composite materials depend on the material being
    % deleted
    comp_depends = false(1,numel(app.matdb));  % remember material dependencies
    comp_log = app.composite & (~k);  % logical index of composite materials
    comp_num = find(comp_log==1);  % numeric index of composite materials
    for c = comp_num  % for each composite material
        mat = app.matdb(c);  % get its data
        layernames = cell(1,mat.uniqueLayers);  % get the layer names
        for j = 1:mat.uniqueLayers
            layernames{j} = mat.layer(j).layerName;
        end
        % determine if this composite uses the material being deleted
        UsesMat = strcmp(matname,layernames);
        if any(UsesMat)
            comp_depends(c) = true;
        end
    end
    %transpose(app.matlist(comp_depends))
    
    % jcb: ToDo - also need to determine how material is used in model
    
    % confirm delete with user
    if any(comp_depends)
        msg=sprintf('Delete material "%s" and the following dependent materials?',matname);
        for j = find(comp_depends==1)
            msg=sprintf('%s\n    "%s"',msg,app.matlist{j});
        end
        ButtonName = questdlg(msg,'Confirm Delete','Delete','Cancel','Cancel');
    else
        msg=sprintf('Delete material "%s"?',matname);
        ButtonName = questdlg(msg,'Confirm Delete','Delete','Cancel','Cancel');
    end
    
    % proceed with delete if user confirms
    if strcmp(ButtonName,'Delete')
            toDelete = k | comp_depends;
            app.matdb(toDelete) = [];  % delete materials
            
            % write out the changes to the material database file
            writeMatDB(app.matdb,app.matdb_path);
            refresh_matdb(app)
            
            guidata(app.fh,app);  % store the application data as guidata
            if app.debugging
                assignin('base','app_mat',app);
            end
            
            refresh_matlistboxes(app.fh);
    end
        
end

function cb_matsave(cbo,~)
% this function is called by another script to save changes to
% a material definition
    subapp = guidata(cbo);
    app = guidata(subapp.caller); 
    
    % find the material in the material list
    k = strcmp(subapp.matname,app.matlist);
    % check if the new name is already in the material list
    k2 = strcmp(subapp.mat.name,app.matlist);
    if any(k2) && any(k2~=k)
        % if the new name was found and it's not the original name
        warndlg('That material name already exists.','Operation Not Permitted')
        return;
    end
    if isempty(subapp.mat.name)
        warndlg('The material name cannot be empty.','Operation Not Permitted')
        return;
    end
    if ~any(k)
        % if the original name was not found, it should be a new material
        k = numel(app.matdb) + 1;  % append to end
        % note: here k changes from logical indexing to numerical indexing
    end
    % save the modified material data into the database
    app.matdb(k) = subapp.mat;
    % propagate changes to other parts of the application
    % jcb:  ToDo name changes affect composite materials
    % jcb:  ToDo material changes affect stored data in blade model
    
    % write out the changes to the material database file
    writeMatDB(app.matdb,app.matdb_path);
    refresh_matdb(app);
    
    guidata(app.fh,app);  % store the application data as guidata
    if app.debugging
        assignin('base','app_mat',app);
    end
    
    refresh_matlistboxes(app.fh);
end


function cb_close(cbo,~)
    app = guidata(cbo);  % load application data
    position = get(app.fh,'Position');
    settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
    if ~isequal(settings.xy_materials,position(1:2))
        settings.xy_materials = position(1:2);
        writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
    end
    closereq;
end

function refresh_matdb(app)
    callapp = guidata(app.caller);

    callapp.matdb = app.matdb;
    callapp.complist = {};
    for k=1:numel(app.matdb)
        if isequal(app.matdb(k).type,'composite')
            callapp.complist{end+1} = app.matdb(k).name;
        end
    end
    guidata(app.caller,callapp);
    
    % populate the list of composite materials
    h = findobj(app.caller,'Tag','DPmaterial');
    val = get(h,'Value');
    str = get(h,'String');
    n = strcmp(callapp.complist,str{val});
    if any(n)
        set(h,'String',callapp.complist,'Value',find(n==1));
    else
        set(h,'String',callapp.complist,'Value',1);
    end 
    h = findobj(app.caller,'Tag','SWmaterial');
    val = get(h,'Value');
    str = get(h,'String');
    n = strcmp(callapp.complist,str{val});
    if any(n)
        set(h,'String',callapp.complist,'Value',find(n==1));
    else
        set(h,'String',callapp.complist,'Value',1);
    end 
end