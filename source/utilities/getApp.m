function app=getApp(varargin)
%NUMAD_MAIN  Main user interface for NuMAD
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_main constructs the main user interface for NuMAD.
%   
%   NuMAD_main(FILENAME) opens a project for editing
%      FILENAME  = filename of a NuMAD project file
%

%===== INITIALIZATION =====================================================

USING_MATLAB_COMPILER = false;
if nargin > 1
    app.batchrun = true;
else
    app.batchrun = false;
end
app.debugging = 0;  % with debugging==1, the data structure is sent to the 
                    % workspace after the gui is constructed
app.ReleaseVersion = 'v3.0.0';
if USING_MATLAB_COMPILER
    app.numadpath = '.';  % point to current directory
else
    numadmain = which('NuMAD_main'); % THIS CAUSES PROBLEMS IN COMPILED VERSION
    if isempty(numadmain)
        errordlg('Could not discover NuMAD path.','error');
        error('Could not discover NuMAD path.');
    end
    [app.numadpath,~,~] = fileparts(numadmain);
end

% ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% locate the user's appdata (windows) or home directory (unix)
% and check for the NuMAD folder
% pc: %APPDATA%\NuMAD\
% unix: $HOME/.numad/
% mac:   ???
parFlag = false; 

global ansysPath
global bmodesPath
global precompPath

parID = gcp('nocreate'); workID = getCurrentWorker;
if ~isempty(parID) || ~isempty(workID)% running in parallel, 
    parFlag = true;
%% EMA Original:
%         
% elseif ispc
%% Changed to:
end

if ispc
%% End
    settingsFile = which('numad');
    [baseFolder,~,~] = fileparts(settingsFile);
    app.userpath = fullfile(baseFolder, 'settings');
    if ~isdir(app.userpath)
        mkdir(app.userpath)
    end
    numadfolder = 'NuMAD';

elseif isunix
    % get the unix/linux $HOME path
    app.userpath = getenv('HOME');
    numadfolder = '.numad';
    if isequal(0,exist(fullfile(app.userpath,numadfolder),'dir'))
        % NuMAD folder not found in %APPDATA%, try creating directory
        [success,message,messageid] = mkdir(app.userpath,numadfolder);
        if isequal(0,success)
            errordlg('Could not create .numad folder in $HOME','Error');
            error(messageid,message);
            return;
        end
    end 
    app.userpath = fullfile(app.userpath,numadfolder);
    
elseif ismac
    errordlg('Mac is not currently supported by NuMAD.','Error');
    error('Mac is not currently supported by NuMAD.');
    return;
    
else
    errordlg('Your system is not supported by NuMAD.','Error');
    error('Your system is not supported by NuMAD.');
    return;
end

% ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if isequal(2,exist(fullfile(app.userpath,'settings.txt'),'file'))
    % settings file found, read settings
    app.settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
else
    % settings file not found; create default settings
    app.settings.ansys_path = ansysPath;
    app.settings.ansys_product = 'ANSYS';
    app.settings.bmodes_path = bmodesPath;
    app.settings.precomp_path = precompPath;
    app.settings.monitors = '';
    app.settings.xy_main = [20 50 1000 600];
    app.settings.xy_materials = [20 50];
    app.settings.xy_isotropic = [20 50];
    app.settings.xy_orthotropic = [20 50];
    app.settings.xy_composite = [20 50];
    app.settings.xy_activelist = [20 50];
    app.settings.xy_ansys = [20 50];
    app.settings.xy_genline = [20 50];
    app.settings.xy_plot3d = [20 50];
    app.settings.xy_pathconfig = [20 50];
    app.settings.xy_flutterinput = [20 50]; 
    app.settings.xy_bpesegments = [20 50 800 300];
    if ispc
        app.settings.job_path = getenv('USERPROFILE');
    elseif isunix
        app.settings.job_path = getenv('HOME');
    else
        app.settings.job_path = '';
    end
    app.settings.job_name = '';
    writeNuMADsettings(app.settings,fullfile(app.userpath,'settings.txt'));
end
% ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


% attempt to read ANSYS product variable file
if isequal(2,exist(fullfile(app.userpath,'ansys-productvars.txt'),'file'))
    % product variable file found
else
    % product variable file not found; regenerate
    filename = fullfile(app.userpath,'ansys-productvars.txt');
    fid = fopen(filename,'wt');
    if (fid == -1)
        error('Could not open file "%s"',filename);
    end
    try
        fprintf(fid,'%% This file lists the ANSYS product variable names available in the NuMAD dropdown lists.\n');
        fprintf(fid,'%% You may manually edit this file to add/remove entries or change the description as you please.\n');
        fprintf(fid,'%% If you delete this file, the default content will regenerate when NuMAD restarts.\n');
        fprintf(fid,'%% See "INSTALLATION AND LICENSING Chapter 6: Product Variable Table" of ANSYS (14.0) help for a list of valid product variable names.\n');
        fprintf(fid,'%% Any line beginning with ''%%'' is ignored.\n');
        fprintf(fid,'%% NAME    DESCRIPTION\n');
        ansys_products = {'ANE3FL','Multiphysics';
            'ANE3FL1','Multiphysics 1';
            'ANE3FL2','Multiphysics 2';
            'ANE3FL3','Multiphysics 3';
            'ANE3FLDS','Multiphysics/LS-DYNA';
            'ANE3FLDP','Multiphysics/LS-DYNA PrepPost';
            'ANSYS','Mechanical';
            'ANE3','Mechanical/Emag';
            'ANFL','Mechanical/FLOTRAN';
            'ANSYSDS','Mechanical/LS-DYNA';
            'ANE3DS','Mechanical/Emag/LS-DYNA';
            'ANFLDS','Mechanical/CFD-Flo/LS-DYNA';
            'ANSYSDP','Mechanical/LS-DYNA PrepPost';
            'ANE3DP','Mechanical/Emag/LS-DYNA PrepPost';
            'ANFLDP','Mechanical/CFD-Flo/LS-DYNA PrepPost';
            'ANCFX','Mechanical/CFD-Flo';
            'STRUCT','Structural';
            'STE3','Structural/Emag';
            'STFL','Structural/FLOTRAN';
            'STE3FL','Structural/Emag/CFD-Flo';
            'STRUCTDS','Structural/LS-DYNA';
            'STE3DS','Structural/Emag/LS-DYNA';
            'STFLDS','Structural/CFD-Flo/LS-DYNA';
            'STE3FLDS','Structural/Emag/CFD-Flo/LS-DYNA';
            'STRUCTDP','Structural/LS-DYNA PrepPost';
            'STE3DP','Structural/Emag/LS-DYNA PrepPost';
            'STFLDP','Structural/CFD-Flo/LS-DYNA PrepPost';
            'STE3FLDP','Structural/Emag/CFD-Flo/LS-DYNA PrepPost';
            'STCFX','Structural/CFD-Flo';
            'PRF','Professional NLT';
            'PRFNLS','Professional NLS';
            'PRFE3','Professional NLT Emag';
            'PRFFL','Professional/FLOTRAN';
            'PRFE3FL','Professional/Emag/FLOTRAN';
            'AA_A','ANSYS Academic Associate';
            'AA_R','ANSYS Academic Research';
            'AA_R_ME','ANSYS Academic Research Mechanical';
            'AA_T_ME','ANSYS Academic Teaching Mechanical'};
        for krow=1:size(ansys_products,1)
            fprintf(fid,'%-8s  %s\n',ansys_products{krow,:});
        end
    catch ME
        fclose(fid);
        rethrow(ME);
    end
    fclose(fid);
end



units = unitconst();
%===== GUI CONSTRUCTION ===================================================
mp = get(0,'MonitorPositions');
if ~strcmp(app.settings.monitors,sprintf('%d ',mp'))
    % monitor settings have changed, reset default window locations
    app.settings.monitors = sprintf('%d ',mp');
    if mp(1,3) >= 1680
        app.settings.xy_main = [20 50 1400 800];
    elseif mp(1,3) >= 1024
        app.settings.xy_main = [20 50 1000 600];
    else
        app.settings.xy_main = [20 50 1000 600];
    end
    app.settings.xy_materials = [20 50];
    app.settings.xy_isotropic = [20 50];
    app.settings.xy_orthotropic = [20 50];
    app.settings.xy_composite = [20 50];
    app.settings.xy_activelist = [20 50];
    app.settings.xy_ansys = [20 50];
    app.settings.xy_genline = [20 50];
    app.settings.xy_plot3d = [20 50];
    app.settings.xy_pathconfig = [20 50];
    app.settings.xy_bpesegments = [20 50 800 300];
    if ~parFlag %ble
        writeNuMADsettings(app.settings,fullfile(app.userpath,'settings.txt'));
    end
end 
c=cgui([],'init');  c.fig.xy = app.settings.xy_main(1:2);
c=cgui(c,'figure','NuMAD',app.settings.xy_main(3:4),'guiNuMAD');
set(c.fig.h,'CloseRequestFcn',@cb_close);
%set(c.fig.h,'toolbar','figure');
set(c.fig.h,'resize','on');
c=cgui(c,'axes','normalized',[0.2 0 0.8 1]);
%view(45,30); %zoom(1.8);
axis([-150 150 -100 100 -10000 10000]);
view(2); %zoom(1.8);
axis manual off;
app.fh = c.fig.h;
app.ax = c.ax.h;
NuMAD_appdata('init','',app.fh);

set( app.fh, 'WindowButtonMotionFcn', '' );
set( app.fh, 'WindowButtonDownFcn', @wf_buttondown );
set( app.fh, 'WindowButtonUpFcn', @wf_buttonup );
set( app.fh, 'WindowKeyPressFcn', @wf_keypress );
set( app.fh, 'WindowKeyReleaseFcn', @wf_keyrelease );
set( app.fh, 'WindowScrollWheelFcn', @wf_scrollwheel );
bview.hgtf = hgtransform('Parent',app.ax);
bview.SVhgtf = hgtransform('Parent',app.ax);
bview.az = -0.8; % azimuth
bview.el = -0.9; % elevation
bview.trans = zeros(3,1);
bview.SVtrans = zeros(3,1);
bview.movement = zeros(3,1);
bview.scale = 1;
bview.SVscale = 1;
bview.lastclick = ''; % SelectType of last mouse click
bview.modifier = '';
bview.function = '';
bview.modifysmd = false;
bview.pointer = getPointerShapes;
NuMAD_appdata('set','bview',bview);

% %jcb: until we create our own toolbar, modify the standard figure toolbar
% htoolbar = findall(app.fh,'Type','uitoolbar');
% hbutton = findall(htoolbar,'Type','uipushtool');
% for k=1:numel(hbutton)
%     switch get(hbutton(k),'TooltipString')
%         case 'New Figure'
%             %jcb: @cb_newmodel does not check if the current model has been
%             %     saved, so for now don't display this button to avoid
%             %     accidental dumping of work
%             delete(hbutton(k));
%             %set(hbutton(k),'ClickedCallback',@cb_newmodel,'TooltipString','New Model')
%         case 'Open File'
%             %delete(hbutton(k));
%             set(hbutton(k),'ClickedCallback',@cb_openmodel,'TooltipString','Open Model...')
%         case 'Save Figure'
%             %delete(hbutton(k));
%             set(hbutton(k),'ClickedCallback',@cb_savemodel,'TooltipString','Save Model')
%         case 'Print Figure'
%             %delete(hbutton(k));
%             set(hbutton(k),'ClickedCallback',@cb_generate,'TooltipString','Generate ANSYS File')
%         case {'Hide Plot Tools',...
%               'Show Plot Tools and Dock Figure'}
%             delete(hbutton(k));
%     end
% end
% hbutton = findall(htoolbar,'Type','uitoggletool');
% for k=1:numel(hbutton)
%     switch get(hbutton(k),'TooltipString')
%         case {'Insert Legend','Insert Colorbar','Link Plot','Edit Plot'}
%             delete(hbutton(k));
%     end
% end
% delete(findall(htoolbar,'Type','uitogglesplittool'));

f=uimenu('label','File');
  uimenu(f,'label','New Model','callback',@cb_newmodel);
  uimenu(f,'label','Open Model...','callback',@cb_openmodel);
  uimenu(f,'label','Save','callback',@cb_savemodel,'separator','on');
  uimenu(f,'label','Save As...','callback',@cb_saveasmodel);
  uimenu(f,'label','Configure Paths...','callback',{@NuMAD_pathconfig,@update_paths},'separator','on');
  uimenu(f,'label','Convert Legacy NuMAD...','callback',@(~,~) NuMAD_converter);
  uimenu(f,'label','XLS-to-NuMAD...','callback',@cb_xls2nmd);
  uimenu(f,'label','Exit','callback',@cb_close,'separator','on');
f=uimenu('label','Blade','tag','uimenu_blade','enable','off');
  uimenu(f,'label','Reference Line','tag','uimenu_genline','enable','off','callback',@NuMAD_genline);
  g=uimenu(f,'label','Rotation','tag','uimenu_rotation','enable','off');
    uimenu(g,'label','Clockwise','tag','uimenu_cw','callback',@cb_bladerotation);
    uimenu(g,'label','Counter-Clockwise','tag','uimenu_ccw','callback',@cb_bladerotation);
f=uimenu('label','View','tag','uimenu_view','enable','off');
  uimenu(f,'label','Translate tool [T]','tag','uimenu_viewtool','callback',@cb_viewtool,'checked','on');
  uimenu(f,'label','Rotate 3D tool [R]','tag','uimenu_viewtool','callback',@cb_viewtool);
  uimenu(f,'label','Zoom tool [Z]','tag','uimenu_viewtool','callback',@cb_viewtool);
  uimenu(f,'label','View Upper (LP)','callback',@cb_orientation,'separator','on');
  uimenu(f,'label','View Lower (HP)','callback',@cb_orientation);
  uimenu(f,'label','View L.E.','callback',@cb_orientation);
  uimenu(f,'label','View T.E.','callback',@cb_orientation);
  uimenu(f,'label','View Root','callback',@cb_orientation);
  uimenu(f,'label','View Tip','callback',@cb_orientation);
  uimenu(f,'label','Reset View','callback',@cb_orientation,'separator','on');
f=uimenu('label','Materials');
  uimenu(f,'label','Material Database','callback',@NuMAD_materials);
  uimenu(f,'label','Material Palette','callback',{@NuMAD_activelist,@update_activelist});
f=uimenu('label','ANSYS','tag','uimenu_ansys','enable','off');
  uimenu(f,'label','Output Options','callback',{@NuMAD_ansys,@update_ansys});
  uimenu(f,'label','Generate Now','callback',@cb_generate);
f=uimenu('label','Plot3D','tag','uimenu_plot3d','enable','off');
  uimenu(f,'label','Output Options','callback',@NuMAD_plot3d);
  uimenu(f,'label','Generate .p3d','callback',@cb_gen_plot3d);
f=uimenu('label','Advanced','tag','uimenu_advanced','enable','off');
  g=uimenu(f,'label','BPE');
  uimenu(g,'label','Set BPE segments','callback',@BPESegments);
  uimenu(g,'label','Run BPE-to-FAST','callback',@cb_runbpe);
  g=uimenu(f,'label','PreComp');
  uimenu(g,'label','Run PreComp-to-FAST','callback',@cb_runprecomp);
  uimenu(f,'label','Flutter Tool','callback',@FlutterToolInput);
f=uimenu('label','Help');
  uimenu(f,'label','User Guide','callback',{@cb_help,'userguide'});
  uimenu(f,'label','User Guide (web)','callback',{@cb_help,'userguideweb'});
  uimenu(f,'label','About NuMAD','callback',{@cb_help,'aboutnumad'});

% f=uimenu('label','Window');
%   uimenu(f,'label','Bring all onscreen','callback','','enable','off');

c=cgui(c,'panel','Station Parameters',[10 200 260 410]);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 120     250]);
c=cgui(c,'select','Airfoil:',{''},'AirfoilName',@cb_airfoil);
c=cgui(c,'select','TE:',{'round','sharp','flat'},'TEtype',@cb_TEtype);
c=cgui(c,'colEdges','',[10 120 200 250]);
c=cgui(c,'dnum','Distance from root:',0.0,'LocationZ',@cb_zloc,units.length);
c=cgui(c,'dnum','Chord length:',1.0,'Chord',@cb_dnum,units.length);
c=cgui(c,'dnum','Twist of station:',0.0,'DegreesTwist',@cb_dnum,units.angle);
c=cgui(c,'num','Normalized X offset:',0.3,'Xoffset',@cb_num);
c=cgui(c,'num','Aerodynamic Center:',0.25,'AeroCenter',@cb_num);
c=cgui(c,'colEdges','',[10         250]);
c=cgui(c,'push','Modify Skin Material Divisions','','modify_smd',@cb_modifySMD);
c=cgui(c,'colEdges','',[10   100      ]);
c=cgui(c,'push','<= Prev. Station','','station_previous',@cb_stationprevious);
c=cgui(c,'vpos','previous');
c=cgui(c,'colEdges','',[   100 130 160   ]);
c=cgui(c,'display','Sta:','','station_current','');
c=cgui(c,'vpos','previous');
c=cgui(c,'colEdges','',[     160   250]);
c=cgui(c,'push','Next Station =>','','station_next',@cb_stationnext);
c=cgui(c,'colEdges','',[10         250]);
c=cgui(c,'push','Check Blade Data','','check_model',@cb_checkbladedata);
c=cgui(c,'vpos','next');
c=cgui(c,'push','New Station','','station_new',@cb_stationnew);
c=cgui(c,'push','Done','','station_done',@cb_stationdone);
c=cgui(c,'push','Cancel','','station_cancel',@cb_stationcancel);
c=cgui(c,'vpos','next');
c=cgui(c,'push','Delete Selected Station','','station_delete',@cb_stationdelete);

c=cgui(c,'panel','Skin Material Division Points',[10 10 260 180]);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 60 140]);
c=cgui(c,'select','number:',{''},'DPnumber',@cb_modifySMD_uictrls);
c=cgui(c,'vpos','previous');
c=cgui(c,'colEdges','',[145  235]);
c=cgui(c,'push','Delete DP','','DPbutton',@cb_modifySMD_uictrls);
c=cgui(c,'colEdges','',[10 60 140]);
c=cgui(c,'select','type:',{'single','double','hourglass','flare'},'DPtype',@cb_modifySMD_uictrls);
c=cgui(c,'vpos','previous');
c=cgui(c,'colEdges','',[135 145 235]);
c=cgui(c,'select','',{''},'DPmaterial',@cb_modifySMD_uictrls);
c=cgui(c,'colEdges','',[10 60 190]);
c=cgui(c,'select','surface:',{'Lower (HP)','Upper (LP)'},'DPsurf',@cb_modifySMD_uictrls);
c=cgui(c,'colEdges','',[10 100 190]);
c=cgui(c,'num','% chord:','','DPchord',@cb_modifySMD_uictrls);
c=cgui(c,'num','chordal distance:','','DPcdist',@cb_modifySMD_uictrls);
c=cgui(c,'num','surface distance:','','DPsdist',@cb_modifySMD_uictrls);

c=cgui(c,'panel','Shear Webs',[280 10 240 180]);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 60 140]);
c=cgui(c,'select','number:',{''},'SWnumber',@cb_modifySW_uictrls);
c=cgui(c,'vpos','previous');
c=cgui(c,'colEdges','',[145  235]);
c=cgui(c,'push','Delete SW','','SWbutton',@cb_modifySW_uictrls);
c=cgui(c,'colEdges','',[10 60 235]);
c=cgui(c,'select','material:',{''},'SWmaterial',@cb_modifySW_uictrls);
c=cgui(c,'colEdges','',[10 80 140]);
c=cgui(c,'display','station:','','SWstation1','');
c=cgui(c,'select','Upper DP:',{''},'SWudp1',@cb_modifySW_uictrls);
c=cgui(c,'select','Lower DP:',{''},'SWldp1',@cb_modifySW_uictrls);
c=cgui(c,'vpos','previous');
c=cgui(c,'vpos','previous');
c=cgui(c,'vpos','previous');
c=cgui(c,'colEdges','',[140 150 210]);
c=cgui(c,'display','','','SWstation2','');
c=cgui(c,'select','',{''},'SWudp2',@cb_modifySW_uictrls);
c=cgui(c,'select','',{''},'SWldp2',@cb_modifySW_uictrls);

app.sm_cmenu = uicontextmenu();  % create context menu for surface material modification
app.dp_cmenu = uicontextmenu();  % create context menu for delineation point modification
h=uimenu(app.dp_cmenu,'Label','DP type','Callback','');
    uimenu(h,'Label','Single','Callback',@cb_dpmodify);
    uimenu(h,'Label','Double','Callback',@cb_dpmodify);
    uimenu(h,'Label','Hourglass','Callback',@cb_dpmodify);
    uimenu(h,'Label','Flare','Callback',@cb_dpmodify);
uimenu(app.dp_cmenu,'Label','Move DP','Callback',@cb_dpmove);
uimenu(app.dp_cmenu,'Label','New DP','Callback','','Enable','off');

[app.afdb app.aflist] = readAirfoilDB(fullfile(app.numadpath,'airfoils'));
app.n_panels = 200;
app.afdb = resampleAirfoilDB(app.afdb,app.n_panels,'cosine');
app.matdb = readMatDB(fullfile(app.numadpath,'MatDBsi.txt'));
app.complist = {};
for k=1:numel(app.matdb)
    if isequal(app.matdb(k).type,'composite')
        app.complist{end+1} = app.matdb(k).name;
    end
end

app.active.list = {'**UNSPECIFIED**'};
%app.active.color = {[1 0 0]};
app.selectenabled = false;
app.modifysmd = false;
app.settings.job_name = '';

app.colors = {[255   0   0]/255;  % red
              [  0 102 255]/255;  % blue
              [  0 204   0]/255;  % green
              [204 102   0]/255;  % brown
              [255 255   0]/255;  % yellow
              [  0 153 255]/255;  % sky blue
              [255 153   0]/255;  % orange
              [  0 153   0]/255;  % dark green
              [204 102 255]/255;  % purple
              [204 153   0]/255;  % light brown
              [  0 255 255]/255;  % cyan
              [255  51   0]/255;  % light red
              [  0 255   0]/255;  % bright green
              [204   0   0]/255;  % dark red
              [255 204 102]/255;  % tan
              [  0 102   0]/255;  % forest green
              [255   0 255]/255;  % magenta
              [102 102 102]/255;  % gray - 40%
              [255   0 102]/255;  % name? (pinkish red)
              [102 153   0]/255}; % olive green


handles.func.draw_stations = @draw_stations;
NuMAD_appdata('set','handles',handles);

blade.BladeRotation = 'ccw';
NuMAD_appdata('set','blade',blade);
guidata(app.fh,app);  % store the application data as guidata

set(app.ax,'ButtonDownFcn',@bdf_stationselect);
gui_disable(app.fh,...
    {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
     'Chord','Chord_units','Xoffset','AeroCenter','LocationZ',...
     'LocationZ_units','modify_smd','station_previous','station_next',...
     'station_new','station_done','station_cancel','station_delete','check_model'});
gui_disable(app.fh,{'DPnumber','DPbutton','DPtype','DPmaterial','DPsurf','DPchord','DPcdist','DPsdist'});
gui_disable(app.fh,{'SWnumber','SWbutton','SWmaterial',...
            'SWstation1','SWudp1','SWldp1','SWstation2','SWudp2','SWldp2'});

if app.debugging
    %assignin('base','c',c);
    assignin('base','app',app);
end

% handle input arguments
if nargin > 0
    if exist(varargin{1},'file')
        [pn,fn,ext] = fileparts(varargin{1});
        if isempty(pn)
            pn = pwd;
        end
        openmodel(pn,[fn,ext],app.fh);
    else
        warndlg(sprintf('Could not find input file: "%s"',varargin{1}),'File Not Found');
    end
end
if nargin > 1
    set(app.fh,'Visible','off');
    try
    batchrun = lower(varargin{2});
    switch batchrun
        % second argument specifies type of batch run
        case {'generate','shell7','ansys'}
            app = guidata(app.fh);  % load application data
            if isequal(batchrun,'shell7');
                app.ansys.dbgen = 0;
                app.ansys.shell7gen = 1;
            end
            if isequal(batchrun,'ansys');
                app.ansys.dbgen = 1;
                app.ansys.shell7gen = 1;
            end
            if nargin > 2
                % third argument is elementsize
                app.ansys.elementsize = varargin{3};
            end
            guidata(app.fh,app);  % store the application data as guidata
            cb_generate(app.fh,[]);
        case 'precomp'
            app = guidata(app.fh);  % load application data
            if nargin > 2
                app.batchmodelist = varargin{3};
            else
                error('Batch run of PreComp-to-FAST requires mode list as third argument.');
            end
            guidata(app.fh,app);  % store the application data as guidata
            cb_runprecomp(app.fh,[])
    end
    close(app.fh); % close NuMAD after finishing batch runs
    catch ME
        set(app.fh,'Visible','on');
        rethrow(ME)
    end
end

end %END GUI CONSTRUCTION (main function)


%===== UTILITY & CALLBACK FUNCTIONS =======================================
function cellIndex = findCellIndex(cellArray,cellValue)
    cellIndex = find(cellfun(@(c) c==cellValue,cellArray));
end

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

function gui_alter(fh,param,setting,tags)
    for k=1:numel(tags)
        h = findobj(fh,'tag',tags{k});
        set(h,param,setting);
    end
end

function uictrl = get_uictrl(fh,tag)
    uictrl.h = findobj(fh,'Tag',tag);
    uictrl.string = get(uictrl.h,'String');
    uictrl.value = get(uictrl.h,'Value');
    uictrl.userdata = get(uictrl.h,'UserData');
    if strcmp(get(uictrl.h,'Style'),'popupmenu')
        uictrl.select = uictrl.string{uictrl.value};
    end
end

function cb_savemodel(cbo,~)
    app = guidata(cbo);  % load application data
    app.blade = NuMAD_appdata('get','blade');

    if isempty(app.settings.job_name)
        % job_name not set, do saveas
        cb_saveasmodel(cbo,[]);
    else
        filename = fullfile(app.settings.job_path,app.settings.job_name);
        writeNuMADinput(app,filename);
        % notify user that model was successfully saved
        msgbox(sprintf('Model saved to %s\n',filename),'Save successful');
    end
end

function cb_saveasmodel(cbo,~)
    app = guidata(cbo);  % load application data
    app.blade = NuMAD_appdata('get','blade');

%     [fn,pn] = uiputfile( ...
%         {'*.nmd', 'NuMAD blade (*.nmd)';...
%          '*.*',   'All files (*.*)'},...
%          'Save As',app.settings.job_path);

    [fn,pn] = uiputfile( ...
        {'*.*',   'Enter project name'},...
         'Save As - enter project name',app.settings.job_path);
    if isequal(fn,0)
        % the user canceled file selection
        % jcb: should we notify user that model was not saved?
        return
    end
    
    % check if directory already exists
    if exist(fullfile(pn,fn),'dir')
        warndlg('Project name already exists.  Please choose another.')
        cb_saveasmodel(cbo,[]);
        return;  % make sure we end here after second call finishes
    end
    % create project directory and report any errors
    [success,message,messageid] = mkdir(pn,fn);
    if ~success
        errordlg(message,'Error');
        return;
    end
    % assemble the filename and write the NuMAD file
    pn = fullfile(pn,fn);  % add the project name to the path
    fn = [fn '.nmd'];  % add the NuMAD extension to the file    
    filename = fullfile(pn,fn);
    writeNuMADinput(app,filename);
    % prepare to copy the material and airfoil databases
    if isempty(app.settings.job_name)
        parent_pn = app.numadpath;  % copy from NuMAD installation directory
    else
        parent_pn = app.settings.job_path;  % copy from current job
    end
    % copy the material database
    [success,message,messageid] = copyfile(fullfile(parent_pn,'MatDBsi.txt'),pn);
    if ~success
        errordlg(message,'Error');
        return;
    end
    % copy the airfoil database
    [success,message,messageid] = copyfile(fullfile(parent_pn,'airfoils'),fullfile(pn,'airfoils'));
    if ~success
        errordlg(message,'Error');
        return;
    end
    % notify user that model was successfully saved
    msgbox(sprintf('Model saved to %s\n',filename),'Save successful');
    
    % start working with the new filename
    app.settings.job_path = fullfile(pn,filesep);
    app.settings.job_name = fn;
    writeNuMADsettings(app.settings,fullfile(app.userpath,'settings.txt'));
    set(app.fh,'Name',sprintf('NuMAD - %s',fullfile(app.settings.job_path,app.settings.job_name)));
    % store application data
    guidata(app.fh,app);  
    if app.debugging
      assignin('base','app',app);
    end
end

function cb_newmodel(cbo,~)
    app = guidata(cbo);  % load application data

%     [fn,pn] = uiputfile( ...
%         {'*.nmd', 'NuMAD blade (*.nmd)'; ...
%          '*.*',   'All files (*.*)'}, ...
%          'New Model');
% 
%     if isequal(fn,0)
%         % the user canceled file selection
%         % jcb: should notify user that model was not saved
%         return
%     end
%     [~,~,ext] = fileparts(fn);
%     if isempty(ext)
%         % add the NuMAD extension if it is missing
%         fn = [fn '.nmd'];
%     end
%     filename = fullfile(pn,fn);
    
    
    % clean up any graphics left over from last model
    delete(get(app.ax,'Children'));
    if isfield(app,'station');
        app = rmfield(app,'station');
    end
    app.hg = {};
    app.swhg = {};
    app.smdp.hg_dp = {};
    app.smdp.hg_sm = {};
    app.smdp.hg_sw = {};
    app.settings.job_name = '';
    set(app.fh,'Name',sprintf('NuMAD - %s','new model'));
    app.active.list = {'**UNSPECIFIED**'};
    %app.active.color = {[1 0 0]};
    
    % load the master airfoil database
    [app.afdb app.aflist] = readAirfoilDB(fullfile(app.numadpath,'airfoils'));
    app.afdb = resampleAirfoilDB(app.afdb,app.n_panels,'cosine');
    
    app.matdb = readMatDB(fullfile(app.numadpath,'MatDBsi.txt'));  % use master database
    app.complist = {'**UNSPECIFIED**'};
    for k=1:numel(app.matdb)
        if isequal(app.matdb(k).type,'composite')
            app.complist{end+1} = app.matdb(k).name;
        end
    end
    
    % create some initial stations for new model
    station.AirfoilName = 'circular';
    station.TEtype = 'round';
    station.DegreesTwist = 0;
    station.LocationZ = 0;
    station.Xoffset = 0.3;
    station.AeroCenter = 0.25;
    station.Chord = 1.0;
    n = strcmp(station.AirfoilName,app.aflist);
    station.coords = app.afdb(n).coords;
    station.dp = transpose([-1.0 -0.3 0.3 1.0]);
    station.dptype = {'single';
                      'single';
                      'single';
                      'single'};
    station.dpmaterial = {'**UNSPECIFIED**';
                          '**UNSPECIFIED**';
                          '**UNSPECIFIED**';
                          '**UNSPECIFIED**'};
    station.sm = {'**UNSPECIFIED**';
                  '**UNSPECIFIED**';
                  '**UNSPECIFIED**'};
    station.xy = app.afdb(n).xy;
    station.LE = app.afdb(n).LE;
    station.s = app.afdb(n).s;
    station.c = app.afdb(n).c;
    
    app.station(1) = station;
    station.LocationZ = 10.0;
    app.station(2) = station;
    
    %app.shearweb = struct('Material',[],'BeginStation',[],'EndStation',[],'Corner',[]);
    app.shearweb = struct([]);
    
    app.ansys = struct('BoundaryCondition', 'cantilevered',...
        'ElementSystem', '281',...
        'MultipleLayerBehavior', 'distinct',...
        'meshing', 'elementsize',...
        'smartsize', 5,...
        'elementsize', 0.0500,...
        'shell7gen', 1,...
        'dbgen', 0);
    
    fcopts = {'EMAX','SMAX','TWSI','TWSR','HFIB','HMAT','PFIB','PMAT',...
    'L3FB','L3MT','L4FB','L4MT','USR1','USR2','USR3','USR4',...
    'USR5','USR6','USR7','USR8','USR9'};
    app.ansys.FailureCriteria = cell(numel(fcopts),2);
    app.ansys.FailureCriteria(:,1) = fcopts';
    app.ansys.FailureCriteria(:,2) = deal({false});
    
    
    app.BladeRotation = 'ccw';
    blade.BladeRotation = app.BladeRotation;
    blade.PresweepRef = struct('method','normal','table',[0 0 0],'pptype','poly');
    blade.PrecurveRef = struct('method','shear','table',[0 0 0],'pptype','poly');
    
    % store application data
    guidata(app.fh,app);
    NuMAD_appdata('set','blade',blade);
    if app.debugging
      assignin('base','app',app);
    end
    
    loadmodel(cbo,[]);  % finish loading required data and display model
end

function cb_openmodel(cbo,~)
    app = guidata(cbo);  % load application data
    if isequal(0,exist(app.settings.job_path,'dir'))
        % job_path not valid, use user's HOME
        if ispc
            app.settings.job_path = getenv('USERPROFILE');
        elseif isunix
            app.settings.job_path = getenv('HOME');
        else
            app.settings.job_path = '';
        end
    end
        
    [fn,pn] = uigetfile( ...
        {'*.nmd', 'NuMAD blade (*.nmd)';...
         '*.*',   'All files (*.*)'},...
         'Open a NuMAD blade model',app.settings.job_path);

    if isequal(fn,0)
        % the user canceled file selection
        return
    else
        openmodel(pn,fn,cbo);
    end
    
end

function openmodel(pn,fn,cbo)
    % load application data
    app = guidata(cbo);  
    app.settings.job_path = fullfile(pn,filesep);
    app.settings.job_name = fn;
    writeNuMADsettings(app.settings,fullfile(app.userpath,'settings.txt'));
    set(app.fh,'Name',sprintf('NuMAD - %s',fullfile(app.settings.job_path,app.settings.job_name)));
    % clean up any graphics left over from last model
    delete(get(app.ax,'Children'));
    app.hg = {};
    app.swhg = {};
    app.smdp.hg_dp = {};
    app.smdp.hg_sm = {};
    app.smdp.hg_sw = {};
    input_file = fullfile(pn,fn);
	[app.station app.shearweb app.active app.ansys app.BladeRotation blade app.plot3d app.flutterInput] = readNuMADinput(input_file);
    app.station = resampleAirfoilDB(app.station,app.n_panels,'cosine');
    [app.afdb app.aflist] = readAirfoilDB(fullfile(pn,'airfoils'));  % use local database
    app.afdb = resampleAirfoilDB(app.afdb,app.n_panels,'cosine');
    
    kmiss=1;
    for k=1:length(app.station)
        if ~any(strcmp(app.station(k).AirfoilName,app.aflist))
            % airfoil was not found in database
            afdb_missing(kmiss).name = app.station(k).AirfoilName; %#ok<*AGROW>
            afdb_missing(kmiss).reference = sprintf('Airfoil "%s" recreated on %s for blade project %s',...
                app.station(k).AirfoilName, date, strrep(app.settings.job_name,'.nmd',''));
            afdb_missing(kmiss).coords = app.station(k).coords;
            afdb_missing(kmiss).xy = app.station(k).xy;
            afdb_missing(kmiss).TEtype = app.station(k).TEtype;
            afdb_missing(kmiss).LE = app.station(k).LE;
            afdb_missing(kmiss).s = app.station(k).s;
            afdb_missing(kmiss).c = app.station(k).c;
            kmiss=kmiss+1;
            fprintf('"%s" not found in airfoil database.\n',app.station(k).AirfoilName);
        end
    end
    
    if kmiss > 1
        % look for unique names (assuming data is identical)
        [aflist_missing,ia] = unique({afdb_missing.name});
        afdb_missing = afdb_missing(ia);
        number_missing = length(afdb_missing);
        
        % append them to the airfoil database
        app.afdb(end+1:end+number_missing) = afdb_missing;
        app.aflist(end+1:end+number_missing) = aflist_missing;
        
        % ask user what to do about missing airfoils
        %% EMA original:
%         response = questdlg( ...
%             'Some airfoils not found in database (see output/command window for list). Would you like to recreate the files?', ...
%             'Missing Airfoils', ...
%             'Yes','No','Yes');
        %% changed to:
        warning('Some airfoils not found in database');
        response = 'No';
        %%
        if isequal('Yes',response)
            for k=1:number_missing
                afname = sprintf('%s.txt',afdb_missing(k).name);
                fullafname = fullfile(app.settings.job_path,'airfoils',afname);
                try
                    writeNuMADAirfoil(afdb_missing(k).coords,...
                        afdb_missing(k).reference,fullafname);
                    fprintf('"%s" written\n',afname);
                catch ME
                    % notify the user of the error, then skip to next file 
                    fprintf('Error writing "%s": %s\n',afname,ME.message);
                    continue
                end
            end
        end
    end  % END kmiss > 1
    
    app.matdb = readMatDB(fullfile(pn,'MatDBsi.txt'));  % use local database
    app.complist = {'**UNSPECIFIED**'};
    for k=1:numel(app.matdb)
        if isequal(app.matdb(k).type,'composite')
            app.complist{end+1} = app.matdb(k).name;
        end
    end
    
    % store application data
    guidata(app.fh,app);  
    NuMAD_appdata('set','blade',blade);
    if app.debugging
      assignin('base','app',app);
    end
   
    loadmodel(cbo,[]);  % finish loading required data and display model
end
    
function loadmodel(cbo,~)
% This function finishes loading required data and then displays the model.
% Call after cb_newmodel() or cb_openmodel().
    % load application data
    app = guidata(cbo);  
    bview = NuMAD_appdata('get','bview');
    
    blade = NuMAD_appdata('get','blade');
    blade = calcGenLinePP(blade); % update the piecewise polynomial
    NuMAD_appdata('set','blade',blade);
    
    delete(get(app.sm_cmenu,'Children'));
    for k = 1:numel(app.active.list)
        uimenu(app.sm_cmenu,'Label',app.active.list{k},...
                            'Callback',@cb_smmodify);
    end
    
    % populate the list of airfoils
    h = findobj(app.fh,'Tag','AirfoilName');
    set(h,'String',app.aflist,'Value',1);

    % populate the list of composite materials
    h = findobj(app.fh,'Tag','DPmaterial');
    set(h,'String',app.complist,'Value',1);
    h = findobj(app.fh,'Tag','SWmaterial');
    set(h,'String',app.complist,'Value',1);
    
    % create the hgtransform object used for orbit/pan/zoom
    bview.hgtf = hgtransform('Parent',app.ax);
    bview.SVhgtf = hgtransform('Parent',app.ax);
    bview.az = -0.8; % azimuth
    bview.el = -0.9; % elevation
    bview.bladetip = abs(app.station(end).LocationZ);
    bview.scale = 250/(1 + bview.bladetip);
    bview.SVscale = bview.scale;
    midLocZ = 0.5*bview.bladetip;
    bview.trans = [-midLocZ*bview.scale; 0; 0];
    bview.SVtrans = zeros(3,1);
    bview.movement = zeros(3,1);
    Raz = makehgtform('zrotate',bview.az);
    Rel = makehgtform('xrotate',bview.el);
    T = makehgtform('translate',bview.trans(1),bview.trans(2),bview.trans(3));
    S = makehgtform('scale',bview.scale);
    set(bview.hgtf,'Matrix',Rel*Raz*T*S);
    T = makehgtform('translate',bview.SVtrans(1),bview.SVtrans(2),bview.SVtrans(3));
    S = makehgtform('scale',bview.SVscale);
    set(bview.SVhgtf,'Matrix',Rel*Raz*T*S);
    bview.lastclick = ''; % SelectType of last mouse click
    bview.modifier = '';
    bview.function = '';
    bview.modifysmd = false;
    
    % create new graphics handle for generating line
    app.hgGL = line(0,0,'Parent',bview.hgtf,'Color',[.7 .7 .7],'LineStyle',':');
    
    % create new graphics handle for each station
    for k = 1:numel(app.station)
        app.hg{k} = line(0,0,'Parent',bview.hgtf,...
            'Color',[0 0 0],'LineWidth',2.0);
        set(app.hg{k},'ButtonDownFcn',@bdf_stationselect);
    end
    
    % create new graphics handles for each shear web
    for k = 1:numel(app.shearweb)
        app.swhg{k} = line(0,0,0,'Parent',bview.hgtf,...
            'Color',[0 0 1],'LineWidth',1.0);
        app.smdp.hg_sw{k} = line(0,0,0,'Parent',bview.SVhgtf,...
            'Color',[0 0 1],'LineWidth',1.0,'Visible','off');
    end
    
%     % vectorize some station parameters
%     for k = 1:numel(app.station)
%         blade.LocationZ(k)    = app.station(k).LocationZ;
%         blade.Chord(k)        = app.station(k).Chord;
%         blade.Xoffset(k)      = app.station(k).Xoffset;
%         blade.DegreesTwist(k) = app.station(k).DegreesTwist;
%     end
    % sort the stations by LocationZ
    [app,sortorder] = sort_stations(app);
%     blade.LocationZ    = blade.LocationZ(sortorder);
%     blade.Chord        = blade.Chord(sortorder);
%     blade.Xoffset      = blade.Xoffset(sortorder);
%     blade.DegreesTwist = blade.DegreesTwist(sortorder);
    
    % create new graphics handles for surface material segments
    % and delineation points
    for k = 1:3  % three stations displayed when modifying materials
        for j = 1:25  % up to 25 surface divisions allowed
            app.smdp.hg_sm{k,j} = line(0,0,'Parent',bview.SVhgtf,...
                'LineWidth',2);
            app.smdp.hg_dp{k,j} = line(0,0,'Parent',bview.SVhgtf);
        end
    end
    cellfun(@(c) set(c,'Visible','off'), app.smdp.hg_sm);
    cellfun(@(c) set(c,'Marker','o','MarkerSize',6,...
        'MarkerEdgeColor','k','MarkerFaceColor','k','Visible','off'), ...
        app.smdp.hg_dp);
 
    % create new graphics handles for skin areas
    for k = 1:2000  %jcb: need to remove hard limit on number of areas
        app.sahg{k} = surface(zeros(2),zeros(2),zeros(2),'Parent',bview.hgtf,...
            'EdgeColor','none','Visible','off','FaceAlpha',0.8,...
            'ButtonDownFcn',@bdf_skinarea);
    end
    
    % create new graphics handles for skin areas in 3-station view
    for k = 1:100  %jcb: need to remove hard limit on number of areas
        app.SVsahg{k} = surface(zeros(2),zeros(2),zeros(2),'Parent',bview.SVhgtf,...
            'EdgeColor','none','Visible','off','FaceAlpha',0.8,...
            'ButtonDownFcn',@bdf_skinarea);
    end
    
    % initialize application flags and button states
    app.selectenabled = true;
    app.modifysmd = false;
    app.checkpassed = false;
    %   ensure the 'modify_smd' button has the proper initial label
    gui_alter(app.fh,'String','Modify Skin Material Divisions',{'modify_smd'});
    
    % store application data
    guidata(app.fh,app);  
    NuMAD_appdata('set','blade',blade);
    NuMAD_appdata('set','bview',bview);
    if app.debugging
      assignin('base','app',app);
    end
        
%     if strcmp(app.BladeRotation,'cw')
%         view(app.ax,135,30);
%     else
%         view(app.ax,45,30);
%     end
    
    % draw the stations
    draw_stations(app);
    
    % modify the Blade Rotation menu
    if isequal(app.BladeRotation,'cw')
        gui_alter(app.fh,'checked','on',{'uimenu_cw'});
        gui_alter(app.fh,'checked','off',{'uimenu_ccw'});
    else
        gui_alter(app.fh,'checked','off',{'uimenu_cw'});
        gui_alter(app.fh,'checked','on',{'uimenu_ccw'});
    end
    
    % disable most gui controls
    gui_disable(app.fh,...
        {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
        'Chord','Chord_units','Xoffset','AeroCenter',...
        'LocationZ','LocationZ_units','modify_smd',...
        'station_delete'});
    gui_disable(app.fh,...
        {'DPnumber','DPbutton','DPtype','DPmaterial','DPsurf','DPchord',...
        'DPcdist','DPsdist'});
    
    % enable specific gui control
    gui_enable(app.fh,...
        {'station_new','station_previous','station_next','check_model',...
         'uimenu_ansys','uimenu_blade','uimenu_genline','uimenu_rotation',...
         'uimenu_view','uimenu_plot3d','uimenu_advanced'});
end

function [app,sortorder] = sort_stations(app)
    nStations = numel(app.station);
    LocationZ = zeros(1,nStations);
    for k = 1:nStations
        LocationZ(k) = app.station(k).LocationZ;
    end
    [~, sortorder] = sort(LocationZ);
    app.station = app.station(sortorder);
    app.hg = app.hg(sortorder);
    
    deletelist = [];
    for ksw = 1:numel(app.shearweb)
        sw = app.shearweb(ksw);
        sw.BeginStation = find(sortorder==sw.BeginStation);
        sw.EndStation = find(sortorder==sw.EndStation);
        if ((sw.EndStation-sw.BeginStation) ~= 1)
            deletelist(end+1) = ksw;
        else
            app.shearweb(ksw) = sw;
        end
    end
    if ~isempty(deletelist)
        app.shearweb(deletelist) = [];
        delete(app.swhg{deletelist});
        delete(app.smdp.hg_sw{deletelist});
        app.swhg(deletelist) = [];
        app.smdp.hg_sw(deletelist) = [];
    end
end

function draw_stations(app)
    if strcmp(app.BladeRotation,'cw')
        twistFlag = -1;  % cw rotor rotation
    else
        twistFlag = 1;  % ccw rotor rotation
    end
    % draw generating line
    blade = NuMAD_appdata('get','blade');
    LocationZ = sort([app.station(:).LocationZ]);
    LocZ_line = linspace(LocationZ(1), LocationZ(end), 100);
    PresweepRef = ppval(blade.PresweepRef.pp,LocZ_line);
    PrecurveRef = ppval(blade.PrecurveRef.pp,LocZ_line);
    set(app.hgGL,'XData',LocZ_line,...
                 'YData',twistFlag*PresweepRef,...
                 'ZData',          PrecurveRef);
    
    % draw stations
    for k = 1:numel(app.station)
        sta = app.station(k);
        x = (sta.xy(:,1) - sta.Xoffset) * sta.Chord * twistFlag;
        y = (sta.xy(:,2)              ) * sta.Chord;
        twist = twistFlag * sta.DegreesTwist * pi/180;
%         xt{k} = cos(twist) * x - sin(twist) * y;
%         yt{k} = sin(twist) * x + cos(twist) * y;
%         zt{k} = ones(size(x)) * sta.LocationZ;
%         xn{k} = sta.c;
        coords(:,1) = cos(twist) * x - sin(twist) * y;
        coords(:,2) = sin(twist) * x + cos(twist) * y;
        coords(:,3) = zeros(size(x));
        coords(:,4) = ones(size(x));
        
        % use the generating line to translate and rotate the coordinates
        [presweep_rot, precurve_rot] = deal(0);
        if isequal(blade.PresweepRef.method,'normal')
            presweep_slope = ppval(blade.PresweepRef.dpp,sta.LocationZ);
            presweep_rot = atan(presweep_slope*twistFlag);
        end
        if isequal(blade.PrecurveRef.method,'normal')
            precurve_slope = ppval(blade.PrecurveRef.dpp,sta.LocationZ);
            precurve_rot = atan(-precurve_slope);
        end
        transX = twistFlag*ppval(blade.PresweepRef.pp,sta.LocationZ);
        transY = ppval(blade.PrecurveRef.pp,sta.LocationZ);
        R = makehgtform('yrotate',presweep_rot,'xrotate',precurve_rot);
        T = makehgtform('translate',transX,transY,sta.LocationZ);
        coords = coords * R' * T';
        xt{k} = coords(:,1);
        yt{k} = coords(:,2);
        zt{k} = coords(:,3);
        xn{k} = sta.c;
        
        set(app.hg{k},'XData',zt{k},...
                      'YData',xt{k},...
                      'ZData',yt{k});
    end
    
    % draw shear webs
    for ksw = 1:numel(app.shearweb)
        k1 = app.shearweb(ksw).BeginStation;
        c1 = app.shearweb(ksw).Corner(1:2);
        dp1 = app.station(k1).dp(1+c1);
        k2 = app.shearweb(ksw).EndStation;
        c2 = app.shearweb(ksw).Corner(3:4);
        dp2 = app.station(k2).dp(1+c2);
        sw(ksw).xt(1:2) = interp1(xn{k1}(2:end-1),xt{k1}(2:end-1),dp1);
        sw(ksw).xt(3:4) = interp1(xn{k2}(2:end-1),xt{k2}(2:end-1),dp2);
        sw(ksw).yt(1:2) = interp1(xn{k1}(2:end-1),yt{k1}(2:end-1),dp1);
        sw(ksw).yt(3:4) = interp1(xn{k2}(2:end-1),yt{k2}(2:end-1),dp2);
        sw(ksw).zt(1:2) = zt{k1}(1:2);
        sw(ksw).zt(3:4) = zt{k2}(1:2);
        sw(ksw).xt(5) = sw(ksw).xt(1);
        sw(ksw).yt(5) = sw(ksw).yt(1);
        sw(ksw).zt(5) = sw(ksw).zt(1);
        
        set(app.swhg{ksw},'XData',sw(ksw).zt,...
                        'YData',sw(ksw).xt,...
                        'ZData',sw(ksw).yt);
    end
    
    cellfun(@(c) set(c,'Visible','off'), app.sahg);
    cellfun(@(c) set(c,'Visible','off'), app.SVsahg);
end

function draw_materialdivisions(app,activestation)
    if strcmp(app.BladeRotation,'cw')
        twistFlag = -1;  % cw rotor rotation
    else
        twistFlag = 1;  % ccw rotor rotation
    end
    cellfun(@(c) set(c,'UIContextMenu','','Visible','off'), app.smdp.hg_sm);
    cellfun(@(c) set(c,'Visible','off'), app.smdp.hg_dp);
    cellfun(@(c) set(c,'Visible','off'), app.smdp.hg_sw);
    cellfun(@(c) set(c,'Visible','off'), app.sahg);
    cellfun(@(c) set(c,'Visible','off'), app.SVsahg);
    threestations = [-1 0 1] + activestation;
    avgChord = [];
    LocZ = zeros(3,1); zFlag = true(3,1);
    for k = 1:3
        if (threestations(k)<1) || (threestations(k)>numel(app.station))
            zFlag(k) = false;
            continue;  % don't try to plot stations before root or after tip
        end
        avgChord(end+1) = app.station(threestations(k)).Chord;
        LocZ(k) = app.station(threestations(k)).LocationZ;
    end
    avgChord = sum(avgChord)/length(avgChord);
    minLocZ = min(LocZ(zFlag));
    maxLocZ = max(LocZ(zFlag));
    LocZ = (LocZ - LocZ(2))*avgChord/abs(maxLocZ-minLocZ);
    
    for k = 1:3
        if (threestations(k)<1) || (threestations(k)>numel(app.station))
            continue;  % don't try to plot stations before root or after tip
        end
        
        sta = app.station(threestations(k));
        x = (sta.xy(:,1) - sta.Xoffset) * sta.Chord * twistFlag;
        y = (sta.xy(:,2)              ) * sta.Chord;
        twist = twistFlag * sta.DegreesTwist * pi/180;
        xt{k} = cos(twist) * x - sin(twist) * y;
        yt{k} = sin(twist) * x + cos(twist) * y;
        %zt{k} = ones(size(x)) * sta.LocationZ;
        xn{k} = sta.c;
        
        % process the surface materials
        sm_mat = cell(150,1);
        sm_count = 0;
        for ksm=1:length(sta.dp)-1
            if (k==2)
                sm.bdf = @bdf_smselect;
                sm.cmenu = app.sm_cmenu;
                sm.linewidth = 2.0;
                kcolor = rem(ksm-1,numel(app.colors))+1;
                sm.color = app.colors{kcolor};
                sm_mat{ksm} = sta.sm{ksm};
                sm_count = length(sta.dp)-1;
            else
                sm.bdf = '';
                sm.cmenu = '';
                sm.linewidth = 2.0;
                sm.color = [0.7 0.7 0.7];  % display as gray
            end
            if threestations(k)==numel(app.station)
                % if plotting last station ...
                sm.name = '';  % no material is assigned
                sm.color = [0.7 0.7 0.7];  % display as gray
                sm.cmenu = '';  % do not allow right-click
                sm.bdf = '';
            end
            sm.pts = transpose(linspace(sta.dp(ksm),sta.dp(ksm+1),50));
            sm.xt = interp1(xn{k}(2:end-1),xt{k}(2:end-1),sm.pts);
            sm.yt = interp1(xn{k}(2:end-1),yt{k}(2:end-1),sm.pts);
            if isequal(sta.TEtype,'flat') && ksm==1
                % HP side of flatback
                sm.xt = [xt{k}(1:2); sm.xt];
                sm.yt = [yt{k}(1:2); sm.yt];
            elseif isequal(sta.TEtype,'flat') && ksm==length(sta.dp)-1
                % LP side of flatback
                sm.xt = [sm.xt; xt{k}(end-1:end)];
                sm.yt = [sm.yt; yt{k}(end-1:end)];
            end
            %sm.zt = avgChord * ones(size(sm.xt)) * (threestations(k)-activestation);
            sm.zt = LocZ(k) * ones(size(sm.xt));
            set(app.smdp.hg_sm{k,ksm},'XData',sm.zt,...
                'YData',sm.xt,...
                'ZData',sm.yt,...
                'Color',sm.color,...
                'LineWidth',sm.linewidth,...
                'ButtonDownFcn',sm.bdf,...
                'UIContextMenu',sm.cmenu,...
                'Visible','on');
        end
        
        if k==2 && sm_count>0
            legend(app.ax,[app.smdp.hg_sm{k,1:sm_count}],sm_mat(1:sm_count),...
                'Interpreter','none','Visible','on');
        end
        
        % process the delineation points
        DPnumber = get_uictrl(app.fh,'DPnumber');
        dp.xt{k} = interp1(xn{k}(2:end-1),xt{k}(2:end-1),sta.dp);
        dp.yt{k} = interp1(xn{k}(2:end-1),yt{k}(2:end-1),sta.dp);
        %dp.zt{k} = avgChord*(threestations(k)-activestation);
        dp.zt{k} = LocZ(k);
        dp.dptype = sta.dptype;
        for kdp=2:length(sta.dp)-1
            if (k==2) && (kdp==(DPnumber.value+1))
                dp.facecolor = [1 1 0];
            else
                dp.facecolor = [0 0 0];
            end
            switch dp.dptype{kdp}
                case 'single'
                    dp.marker = 'o';
                    dp.markersize = 5;
                case 'double'
                    dp.marker = 'v';
                    dp.markersize = 6;
                case 'hourglass'
                    dp.marker = 'd';
                    dp.markersize = 6;
                case 'flare'
                    dp.marker = '^';
                    dp.markersize = 6;
            end
            set(app.smdp.hg_dp{k,kdp},'XData',dp.zt{k},...
                'YData',dp.xt{k}(kdp),...
                'ZData',dp.yt{k}(kdp),...
                'Marker',dp.marker,...
                'MarkerSize',dp.markersize,...
                'MarkerFaceColor',dp.facecolor,... %'UIContextMenu',app.dp_cmenu,...
                'Visible','on');
        end
    
    end  % end loop over threestations
    
    SWnumber = get_uictrl(app.fh,'SWnumber');
    SWindex = get(SWnumber.h,'UserData');
    for ksw = 1:numel(app.shearweb)
        k1 = app.shearweb(ksw).BeginStation;
        if ~any(k1==threestations(1:2))
            continue;  % skip; this shear web not at these stations 
        end
        c1 = app.shearweb(ksw).Corner(1:2);
        k2 = app.shearweb(ksw).EndStation;
        c2 = app.shearweb(ksw).Corner(3:4);
        sw(ksw).xt(1:2) = dp.xt{k1==threestations}(1+c1);
        sw(ksw).xt(3:4) = dp.xt{k2==threestations}(1+c2);
        sw(ksw).yt(1:2) = dp.yt{k1==threestations}(1+c1);
        sw(ksw).yt(3:4) = dp.yt{k2==threestations}(1+c2);
        sw(ksw).zt(1:2) = dp.zt{k1==threestations};
        sw(ksw).zt(3:4) = dp.zt{k2==threestations};
        sw(ksw).xt(5) = sw(ksw).xt(1);
        sw(ksw).yt(5) = sw(ksw).yt(1);
        sw(ksw).zt(5) = sw(ksw).zt(1);
        
        if (SWnumber.value<=numel(SWindex)) && (ksw==SWindex(SWnumber.value))
            linewidth = 1.0;
            linecolor = [0 0 1];
        else
            linewidth = 1.0;
            linecolor = [0.8 0.8 1];
        end
        set(app.smdp.hg_sw{ksw},'XData',sw(ksw).zt,...
                                'YData',sw(ksw).xt,...
                                'ZData',sw(ksw).yt,...
                                'LineWidth',linewidth,...
                                'Color',linecolor,...
                                'Visible','on');
    end
    
    if isfield(app,'newsw')
        k1 = app.newsw.BeginStation;
        k2 = app.newsw.EndStation;
        c1 = app.newsw.Corner(1:2);
        c2 = app.newsw.Corner(3:4);
        newsw.xt(1:2) = dp.xt{k1==threestations}(1+c1);
        newsw.xt(3:4) = dp.xt{k2==threestations}(1+c2);
        newsw.yt(1:2) = dp.yt{k1==threestations}(1+c1);
        newsw.yt(3:4) = dp.yt{k2==threestations}(1+c2);
        newsw.zt(1:2) = dp.zt{k1==threestations};
        newsw.zt(3:4) = dp.zt{k2==threestations};
        newsw.xt(5) = newsw.xt(1);
        newsw.yt(5) = newsw.yt(1);
        newsw.zt(5) = newsw.zt(1);
        linewidth = 1.0;
        linecolor = [1 0 0];
        set(app.hg_newsw,'XData',newsw.zt,...
            'YData',newsw.xt,...
            'ZData',newsw.yt,...
            'LineWidth',linewidth,...
            'Color',linecolor,...
            'Visible','on');
    end
    
    if isfield(app,'newdp')
        switch app.newdp.type
            case 'single'
                dp.marker = 'o';
                dp.markersize = 5;
            case 'double'
                dp.marker = 'v';
                dp.markersize = 6;
            case 'hourglass'
                dp.marker = 'd';
                dp.markersize = 6;
            case 'flare'
                dp.marker = '^';
                dp.markersize = 6;
        end
        dp.xt = interp1(xn{2}(2:end-1),xt{2}(2:end-1),app.newdp.chord);
        dp.yt = interp1(xn{2}(2:end-1),yt{2}(2:end-1),app.newdp.chord);
        %dp.zt = 0;
        dp.zt = LocZ(2);
        set(app.newdp.hg,'XData',dp.zt,...
                         'YData',dp.xt,...
                         'ZData',dp.yt,...
                         'Marker',dp.marker,...
                         'MarkerSize',dp.markersize,...
                         'MarkerFaceColor',[1 0 0],...
                         'Visible','on');
    end
    
    % install button down functions
%     set(sta.hsm,'ButtonDownFcn',@bdf_selectGraphics);
%     set(sta.hdp,'ButtonDownFcn',@bdf_selectGraphics);
end

function bdf_smselect(cbo,~)
    app = guidata(cbo);   % retrieve application data
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    n = findCellIndex(app.smdp.hg_sm(2,:),cbo);   % get index of surface segment
    
    % checkmark the entry corresponding to the currently selected material
    c = strcmp(app.station(k).sm{n},app.active.list);  % get the index of the material name
    cmenu_items = get(app.sm_cmenu,'Children');
    set(cmenu_items,'Checked','off');
    if any(c)
        set(cmenu_items(flipud(c)),'Checked','on');  
    end
end


function app = draw_skinareas(app)
%     cellfun(@(c) set(c,'UIContextMenu','','Visible','off'), app.smdp.hg_sm);
%     cellfun(@(c) set(c,'Visible','off'), app.smdp.hg_dp);
%     cellfun(@(c) set(c,'Visible','off'), app.smdp.hg_sw);
    cellfun(@(c) set(c,'Visible','off'), app.sahg);
    cellfun(@(c) set(c,'Visible','off'), app.SVsahg);
    blade = NuMAD_appdata('get','blade');
    
    if app.modifysmd
        areaLimit=numel(app.SVsahg);
        active = findobj(gcbf,'Tag','active station');  % get handle of active station
        activestation = findCellIndex(app.hg,active);  % get index of active station
        threestations = [-1 0 1] + activestation;
        % range must be from first station to last station
        outOfRange = (threestations < 1) | (threestations > numel(app.station));
        threestations(outOfRange) = [];
        avgChord = mean([app.station(threestations).Chord]);
        LocZ = [app.station(threestations).LocationZ];
        thisLocZ = app.station(activestation).LocationZ;
        LocZ = (LocZ - thisLocZ)*avgChord/abs(max(LocZ)-min(LocZ));
        stationList = threestations(1:end-1);
    else
        areaLimit=numel(app.sahg);
        stationList = 1:numel(app.SkinAreas);
    end
    
    if strcmp(app.BladeRotation,'cw')
        twistFlag = -1;  % cw rotor rotation
    else
        twistFlag = 1;  % ccw rotor rotation
    end
    areaCounter = 1;
    TOO_MANY_AREAS = false;
    for kStation = stationList
        sta = app.station(kStation);
        thisIsFlatback = isequal(sta.TEtype,'flat');
        x = (sta.xy(:,1) - sta.Xoffset) * sta.Chord * twistFlag;
        y = (sta.xy(:,2)              ) * sta.Chord;
        twist = twistFlag * sta.DegreesTwist * pi/180;
%         ib.xt = cos(twist) * x - sin(twist) * y;
%         ib.yt = sin(twist) * x + cos(twist) * y;
%         ib.zt = ones(size(x)) * sta.LocationZ;
%         ib.xn = sta.c;
        
        % use the generating line to translate and rotate the coordinates
        coords(:,1) = cos(twist) * x - sin(twist) * y;
        coords(:,2) = sin(twist) * x + cos(twist) * y;
        coords(:,3) = zeros(size(x));
        coords(:,4) = ones(size(x));
        [presweep_rot, precurve_rot] = deal(0);
        if isequal(blade.PresweepRef.method,'normal')
            presweep_slope = ppval(blade.PresweepRef.dpp,sta.LocationZ);
            presweep_rot = atan(presweep_slope*twistFlag);
        end
        if isequal(blade.PrecurveRef.method,'normal')
            precurve_slope = ppval(blade.PrecurveRef.dpp,sta.LocationZ);
            precurve_rot = atan(-precurve_slope);
        end
        transX = twistFlag*ppval(blade.PresweepRef.pp,sta.LocationZ);
        transY = ppval(blade.PrecurveRef.pp,sta.LocationZ);
        if app.modifysmd
            R = makehgtform; % identity
            k3 = find(kStation==threestations);
            T = makehgtform('translate',0,0,LocZ(k3));
        else
            R = makehgtform('yrotate',presweep_rot,'xrotate',precurve_rot);
            T = makehgtform('translate',transX,transY,sta.LocationZ);
        end
        coords = coords * R' * T';
        ib.xt = coords(:,1);
        ib.yt = coords(:,2);
        ib.zt = coords(:,3);
        ib.xn = sta.c;
        
        sta = app.station(kStation+1);
        nextIsFlatback = isequal(sta.TEtype,'flat');
        x = (sta.xy(:,1) - sta.Xoffset) * sta.Chord * twistFlag;
        y = (sta.xy(:,2)              ) * sta.Chord;
        twist = twistFlag * sta.DegreesTwist * pi/180;
%         ob.xt = cos(twist) * x - sin(twist) * y;
%         ob.yt = sin(twist) * x + cos(twist) * y;
%         ob.zt = ones(size(x)) * sta.LocationZ;
%         ob.xn = sta.c;

        % use the generating line to translate and rotate the coordinates
        coords(:,1) = cos(twist) * x - sin(twist) * y;
        coords(:,2) = sin(twist) * x + cos(twist) * y;
        coords(:,3) = zeros(size(x));
        coords(:,4) = ones(size(x));
        [presweep_rot, precurve_rot] = deal(0);
        if isequal(blade.PresweepRef.method,'normal')
            presweep_slope = ppval(blade.PresweepRef.dpp,sta.LocationZ);
            presweep_rot = atan(presweep_slope*twistFlag);
        end
        if isequal(blade.PrecurveRef.method,'normal')
            precurve_slope = ppval(blade.PrecurveRef.dpp,sta.LocationZ);
            precurve_rot = atan(-precurve_slope);
        end
        transX = twistFlag*ppval(blade.PresweepRef.pp,sta.LocationZ);
        transY = ppval(blade.PrecurveRef.pp,sta.LocationZ);
        if app.modifysmd
            R = makehgtform; % identity
            k3 = find(kStation==threestations);
            T = makehgtform('translate',0,0,LocZ(1+k3));
        else
            R = makehgtform('yrotate',presweep_rot,'xrotate',precurve_rot);
            T = makehgtform('translate',transX,transY,sta.LocationZ);
        end
        coords = coords * R' * T';
        ob.xt = coords(:,1);
        ob.yt = coords(:,2);
        ob.zt = coords(:,3);
        ob.xn = sta.c;
        
        completeAreas = min([numel(app.SkinAreas(kStation).startIB), numel(app.SkinAreas(kStation).startOB)]);
        for kArea = 1:completeAreas
            sa.pts = linspace(app.station(kStation).dp(app.SkinAreas(kStation).startIB(kArea)),...
                app.station(kStation).dp(app.SkinAreas(kStation).endIB(kArea)),20);
            n = find(strcmp(app.SkinAreas(kStation).Material{kArea},app.matlist)==1);
            %app.SkinAreas(kStation).Material{kArea}
            kcolor = rem(n-1,numel(app.colors))+1;
            sa.color = app.colors{kcolor};
            sa.xt(1,:) = interp1(ib.xn(2:end-1),ib.xt(2:end-1),sa.pts);
            sa.yt(1,:) = interp1(ib.xn(2:end-1),ib.yt(2:end-1),sa.pts);
            sa.zt(1,:) = interp1(ib.xn(2:end-1),ib.zt(2:end-1),sa.pts);
            sa.pts = linspace(app.station(kStation+1).dp(app.SkinAreas(kStation).startOB(kArea)),...
                app.station(kStation+1).dp(app.SkinAreas(kStation).endOB(kArea)),20);
            sa.xt(2,:) = interp1(ob.xn(2:end-1),ob.xt(2:end-1),sa.pts);
            sa.yt(2,:) = interp1(ob.xn(2:end-1),ob.yt(2:end-1),sa.pts);
            sa.zt(2,:) = interp1(ob.xn(2:end-1),ob.zt(2:end-1),sa.pts);
            if (thisIsFlatback && kArea==1)
                sa.xt(1,1) = ib.xt(1);
                sa.yt(1,1) = ib.yt(1);
                sa.zt(1,1) = ib.zt(1);
            end
            if (thisIsFlatback && kArea==numel(app.SkinAreas(kStation).startIB))
                sa.xt(1,end) = ib.xt(end);
                sa.yt(1,end) = ib.yt(end);
                sa.zt(1,end) = ib.zt(end);
            end
            if (nextIsFlatback && kArea==1)
                sa.xt(2,1) = ob.xt(1);
                sa.yt(2,1) = ob.yt(1);
                sa.zt(2,1) = ob.zt(1);
            end
            if (nextIsFlatback && kArea==numel(app.SkinAreas(kStation).startOB))
                sa.xt(2,end) = ob.xt(end);
                sa.yt(2,end) = ob.yt(end);
                sa.zt(2,end) = ob.zt(end);
            end
            %app.sahg{kStation,kArea} = surface(sa.zt,sa.xt,sa.yt,'FaceColor',sa.color,'EdgeColor','none');
            if areaCounter>areaLimit   
                TOO_MANY_AREAS = true;
            else
                if app.modifysmd
                    set(app.SVsahg{areaCounter},...
                        'XData',sa.zt,...
                        'YData',sa.xt,...
                        'ZData',sa.yt,...
                        'FaceColor',sa.color,...
                        'Visible','on');
                else
                    set(app.sahg{areaCounter},...
                        'XData',sa.zt,...
                        'YData',sa.xt,...
                        'ZData',sa.yt,...
                        'FaceColor',sa.color,...
                        'Visible','on');
                    
                end
            end
            areaCounter = areaCounter + 1;
        end
    end
    if TOO_MANY_AREAS
        errordlg('Reached hard-coded limit on number of graphics areas.');
    end
end


% function cb_text(cbo,~)
% % this callback reads text input and updates app data
%     S = guidata(cbo);     % retrieve guidata
%     tag = get(cbo,'tag'); % get the object tag (data struct fieldname)
%     str = get(cbo,'string');  % get the new input
%     S.data.(tag) = str;
%     guidata(S.fh,S);  % store the updated guidata
% end

function cb_num(cbo,~)
% this callback checks numeric input and updates app data
    app = guidata(cbo);     % retrieve application data
    blade = NuMAD_appdata('get','blade');
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    tag = get(cbo,'tag');  % get the object tag (data struct fieldname)
    % get the stored data associated with the fieldname
    val = app.station(k).(tag);  % use dynamic fieldname to get value
    str = get(cbo,'string');     % get the new input
    num = str2double(str);  % try to convert to number
    if isnan(num)
        % input was not valid => revert to stored value
        str = sprintf('%g',val);
        set(cbo,'string',str);
    else
        % input valid => save new value
        app.station(k).(tag) = num;
        %blade.(tag)(k) = app.station(k).(tag);
        guidata(app.fh,app);  % store the updated guidata
        NuMAD_appdata('set','blade',blade);
        if app.debugging
            assignin('base','app',app);
        end
        draw_stations(app);  % redraw the graphics
        if app.modifysmd
            draw_materialdivisions(app,k)
        end
    end
end

function cb_dnum(cbo,~)
% this callback checks dimensional numeric input and updates app data
    app = guidata(cbo);     % retrieve guidata
    blade = NuMAD_appdata('get','blade');
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    style = get(cbo,'style');  % get the uicontrol style
    switch style
        case 'edit'
            % get the object tag (data struct fieldname)
            obj = cbo;
            tag = get(obj,'tag');
            % get the units
            obj_units = findobj(app.fh,'tag',[tag '_units']);
            units = get(obj_units,'UserData');
            currentunit = get(obj,'Value');
            % get the stored data associated with the fieldname
            val = app.station(k).(tag);  % use dynamic fieldname to get value
            val = val * units(currentunit);  % do unit conversion
            str = get(obj,'string');  % get the new input
            num = str2double(str);  % try to convert to number
            if isnan(num)
                % input was not valid => revert to stored value
                str = sprintf('%g',val);
                set(obj,'string',str);
            else
                % input valid => save new value
                app.station(k).(tag) = num / units(currentunit);
                %blade.(tag)(k) = app.station(k).(tag);
                guidata(app.fh,app);  % store the updated guidata
                NuMAD_appdata('set','blade',blade);
                if app.debugging
                    assignin('base','app',app);
                end
                draw_stations(app);  % redraw the graphics
                if app.modifysmd
                    draw_materialdivisions(app,k)
                end
            end
        case 'popupmenu'
            % get the object tag (data struct fieldname)
            obj_units = cbo;
            tag = get(obj_units,'tag');
            tag = tag(1:end-6);  % remove the '_units' ending
            obj = findobj(app.fh,'tag',tag);
            oldunit = get(obj,'Value');
            newunit = get(obj_units,'Value');
            if newunit ~= oldunit
                units = get(obj_units,'UserData');
                val = app.station(k).(tag); % use dynamic fieldname to get value
                val = val * units(newunit);
                str = sprintf('%g',val);
                set(obj,'String',str);
                set(obj,'Value',newunit);
            end
    end
end

function cb_airfoil(cbo,~)
% this callback processes changes to the airfoil dropdown menu
    app = guidata(cbo);   % retrieve application data
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    menulist = get(cbo,'string');
    selection = get(cbo,'value');
    newairfoil = menulist{selection};
    
    n = strcmp(newairfoil,app.aflist);
    app.station(k).coords = app.afdb(n).coords;
    app.station(k).xy = app.afdb(n).xy;
    app.station(k).LE = app.afdb(n).LE;
    app.station(k).c = app.afdb(n).c;
    app.station(k).s = app.afdb(n).s;
    app.station(k).AirfoilName = newairfoil;
    app.station(k).TEtype = app.afdb(n).TEtype;
    
    tag = 'TEtype';
    obj = findobj(app.fh,'Tag',tag);
    list = get(obj,'String');
    n = strcmp(app.station(k).(tag),list);
    set(obj,'Value',find(n==1));
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
    draw_stations(app);  % redraw the graphics
    if app.modifysmd
        draw_materialdivisions(app,k)
    end
end

function cb_zloc(cbo, ~)
    cb_dnum(cbo,[]);
    
    app = guidata(cbo);   % retrieve application data
%    blade = NuMAD_appdata('get','blade');
    [app,sortorder] = sort_stations(app);
%     blade.LocationZ    = blade.LocationZ(sortorder);
%     blade.Chord        = blade.Chord(sortorder);
%     blade.Xoffset      = blade.Xoffset(sortorder);
%     blade.DegreesTwist = blade.DegreesTwist(sortorder);
    
    guidata(app.fh,app);  % store the updated application data
%    NuMAD_appdata('set','blade',blade);
    
    bview = NuMAD_appdata('get','bview');
    bview.bladetip = abs(app.station(end).LocationZ);
    NuMAD_appdata('set','bview',bview);
    
    if app.debugging
        assignin('base','app',app);
    end
end

function cb_TEtype(cbo,~)
% this callback processes changes to the TE type dropdown menu
    app = guidata(cbo);   % retrieve application data
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    
    TEtype = get_uictrl(app.fh,'TEtype');
    app.station(k).TEtype = TEtype.select;
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app_comp',app);
    end
end

function readActiveStation(cbo)
% this function populates the gui controls with the active station's
% parameter values
    % load application data
    app = guidata(cbo);
    % determine which station is active
    k = findCellIndex(app.hg,cbo);
    
    tag = 'AirfoilName';
    obj = findobj(app.fh,'Tag',tag);
    list = get(obj,'String');
    n = strcmp(app.station(k).(tag),list);
    set(obj,'Value',find(n==1));

    tag = 'TEtype';
    obj = findobj(app.fh,'Tag',tag);
    list = get(obj,'String');
    n = strcmp(app.station(k).(tag),list);
    set(obj,'Value',find(n==1));
    
    tag = 'DegreesTwist';
    obj = findobj(app.fh,'Tag',tag);
    obj_units = findobj(app.fh,'Tag',[tag '_units']);
    units = get(obj_units,'UserData');
    currentunit = get(obj,'Value');
    val = app.station(k).(tag); % use dynamic fieldname to get value
    val = val * units(currentunit); 
    str = sprintf('%g',val);
    set(obj,'string',str);
    
    tag = 'Chord';
    obj = findobj(app.fh,'Tag',tag);
    obj_units = findobj(app.fh,'Tag',[tag '_units']);
    units = get(obj_units,'UserData');
    currentunit = get(obj,'Value');
    val = app.station(k).(tag); % use dynamic fieldname to get value
    val = val * units(currentunit); 
    str = sprintf('%g',val);
    set(obj,'string',str);
  
    tag = 'Xoffset';
    obj = findobj(gcbf,'Tag',tag);
    str = sprintf('%g',app.station(k).(tag));
    set(obj,'string',str);
    
    tag = 'AeroCenter';
    obj = findobj(gcbf,'Tag',tag);
    str = sprintf('%g',app.station(k).(tag));
    set(obj,'string',str);

    tag = 'LocationZ';
    obj = findobj(app.fh,'Tag',tag);
    obj_units = findobj(app.fh,'Tag',[tag '_units']);
    units = get(obj_units,'UserData');
    currentunit = get(obj,'Value');
    val = app.station(k).(tag); % use dynamic fieldname to get value
    val = val * units(currentunit); 
    str = sprintf('%g',val);
    set(obj,'string',str);
    
end

function bdf_stationselect(cbo,~)
% this "button down function" responds to mouse clicks on graphics
app = guidata(cbo);  % retrieve the application data
if (app.selectenabled && ~app.modifysmd)
    % look for active (already selected) station
    active = findobj(gcbf,'Tag','active station');
    click = get(gcbf,'SelectionType');  % find out which mouse button was pressed
    switch click
        case {'normal','extend','alt'}  % any mouse button
            if ~isempty(active)
                % station found, deselect the station
                set(active,'Color',[0 0 0],'Tag','');
                gui_alter(gcbf,'String','',{'station_current'});
                gui_disable(gcbf,...
                    {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
                    'Chord','Chord_units','Xoffset','AeroCenter',...
                    'LocationZ','LocationZ_units','modify_smd',...
                    'station_delete'});
            end
            if strcmp(get(cbo,'Type'),'line')
                % if the cbo is a line, make it the active station
                set(cbo,'Color',[0 0.8 0],'Tag','active station');
                newindex = findCellIndex(app.hg,cbo);
                gui_alter(gcbf,'String',num2str(newindex),{'station_current'});
                readActiveStation(cbo);
                gui_enable(gcbf,...
                    {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
                    'Chord','Chord_units','Xoffset','AeroCenter',...
                    'LocationZ','LocationZ_units','modify_smd',...
                    'station_delete'});
            end
            
        case 'open'  % double mouse click (any button)
            if ~isempty(active)
                obj = findobj(app.fh,'Tag','modify_smd');
                cb_modifySMD(obj,[]);
            end
    end
end
end

function bdf_skinarea(cbo,~)
% this "button down function" responds to mouse clicks on graphics
app = guidata(cbo);  % retrieve the application data
click = get(gcbf,'SelectionType');  % find out which mouse button was pressed
if app.modifysmd
    switch click
        case {'alt'}  % right click
            cellfun(@(c) set(c,'Visible','off'), app.SVsahg);
    end
else
    switch click
        case {'alt','open'}  % right click or double click
            cellfun(@(c) set(c,'Visible','off'), app.sahg);
    end
end
end

function cb_stationnew(cbo,~)
% start defining a new station
    app = guidata(cbo);  % retrieve the application data
    bview = NuMAD_appdata('get','bview');
    gui_enable(gcbf,{'station_done','station_cancel'});
    gui_disable(gcbf,{'station_new','station_previous','station_next','check_model'});
    bdf_stationselect(app.ax,[]);   % unselect active station
                
    % create some initial stations for new model
    station.AirfoilName = 'circular';
    station.TEtype = 'round';
    station.DegreesTwist = 0;
    station.LocationZ = 0;
    station.Xoffset = 0.3;
    station.AeroCenter = 0.25;
    station.Chord = 1.0;
    n = strcmp(station.AirfoilName,app.aflist);
    if ~any(n)
        station.AirfoilName = app.afdb(1).name;
        station.TEtype = app.afdb(1).TEtype;
        n(1) = 1;
    end
    station.coords = app.afdb(n).coords;
    station.dp = transpose([-1.0 -0.3 0.3 1.0]);
    station.dptype = {'single';
                      'single';
                      'single';
                      'single'};
    station.dpmaterial = {'**UNSPECIFIED**';
                          '**UNSPECIFIED**';
                          '**UNSPECIFIED**';
                          '**UNSPECIFIED**'};
    station.sm = {'**UNSPECIFIED**';
                  '**UNSPECIFIED**';
                  '**UNSPECIFIED**'};
    station.xy = app.afdb(n).xy;
    station.LE = app.afdb(n).LE;
    station.s = app.afdb(n).s;
    station.c = app.afdb(n).c;
    
    
    new = 1 + numel(app.station);
    app.station(new) = station;
    app.hg{new} = line(0,0,'Parent',bview.hgtf,...
        'Color','red','LineWidth',2.0,'Tag','active station');
    set(app.hg{new},'ButtonDownFcn',@bdf_stationselect);
    guidata(app.fh,app);  % store the updated application data
    
    readActiveStation(app.hg{new});
    app.selectenabled = false;  % disallow selecting other stations
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
    draw_stations(app);  % redraw the graphics
    
    gui_enable(gcbf,...
        {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
        'Chord','Chord_units','Xoffset','AeroCenter',...
        'LocationZ','LocationZ_units'});
    
    gui_alter(gcbf,'BackgroundColor',[1 1 0.8],...
        {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
        'Chord','Chord_units','Xoffset','AeroCenter',...
        'LocationZ','LocationZ_units'});
    
end

function cb_stationdone(cbo,~)
% finish creating a new station
    app = guidata(cbo);  % retrieve the application data
    blade = NuMAD_appdata('get','blade');
    
    gui_disable(gcbf,{'station_done','station_cancel'});
    gui_enable(gcbf,{'station_new','station_delete','modify_smd','station_previous','station_next','check_model'});
    
    gui_alter(gcbf,'BackgroundColor',[1 1 1],...
        {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
        'Chord','Chord_units','Xoffset','AeroCenter',...
        'LocationZ','LocationZ_units'});
    [app,sortorder] = sort_stations(app);
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    set(app.hg{k},'Color',[0 0.8 0]); % change color from red to green
    gui_alter(gcbf,'String',num2str(k),{'station_current'});  % update active station number
%     blade.LocationZ    = blade.LocationZ(sortorder);
%     blade.Chord        = blade.Chord(sortorder);
%     blade.Xoffset      = blade.Xoffset(sortorder);
%     blade.DegreesTwist = blade.DegreesTwist(sortorder);
    
    app.selectenabled = true;  % allow selecting stations again
    guidata(app.fh,app);  % store the updated application data
    NuMAD_appdata('set','blade',blade);
    if app.debugging
        assignin('base','app',app);
    end
    
    bview = NuMAD_appdata('get','bview');
    bview.bladetip = abs(app.station(end).LocationZ);
    NuMAD_appdata('set','bview',bview);
end

function cb_stationcancel(cbo,~)
% cancel creating a new station
    gui_disable(gcbf,{'station_done','station_cancel'});
    gui_enable(gcbf,{'station_new','station_previous','station_next','check_model'});

    gui_alter(gcbf,'BackgroundColor',[1 1 1],...
        {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
        'Chord','Chord_units','Xoffset','AeroCenter',...
        'LocationZ','LocationZ_units'});
    gui_disable(gcbf,...
        {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
        'Chord','Chord_units','Xoffset','AeroCenter',...
        'LocationZ','LocationZ_units'});
    
    cb_stationdelete(cbo,[]);  % delete the new station
end

function cb_stationdelete(cbo,~)
% delete the active (currently selected) station
    app = guidata(cbo);   % retrieve application data
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    bdf_stationselect(app.ax,[]);   % unselect active station
    delete(app.hg{k});  % delete the station's curve
    app.hg(k) = []; % delete the graphic handle 
    app.station(k) = []; % delete the station's data
    % delete associated shear webs
    deletelist = [];
    for ksw = 1:numel(app.shearweb)
        sw = app.shearweb(ksw);
        if (k==sw.BeginStation) || (k==sw.EndStation)
            deletelist(end+1) = ksw;
        elseif (sw.BeginStation > k)
            % for shear webs outboard of station being deleted, shift
            % indices down by one
            app.shearweb(ksw).BeginStation = sw.BeginStation - 1;
            app.shearweb(ksw).EndStation = sw.EndStation - 1;
        end
    end
    if ~isempty(deletelist)
        app.shearweb(deletelist) = [];
        delete(app.swhg{deletelist});
        delete(app.smdp.hg_sw{deletelist});
        app.swhg(deletelist) = [];
        app.smdp.hg_sw(deletelist) = [];
    end
    draw_stations(app);  % redraw the graphics
    
    app.selectenabled = true;  % allow selecting stations again (if new station cancelled)
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
end

function cb_stationprevious(cbo,~)
    app = guidata(cbo);  % retrieve the application data
    if app.selectenabled
        % look for active (already selected) station
        active = findobj(gcbf,'Tag','active station');
        if ~isempty(active)
            % station found, remember station number
            index = findCellIndex(app.hg,active);
            % deselect the station
            set(active,'Color',[0 0 0],'Tag','');
            newindex = index - 1;
            if newindex <= 0
                newindex = numel(app.hg) + newindex;
            end
            newactive = app.hg{newindex};
        else
            newactive = app.hg{end};
            newindex = numel(app.hg);
            gui_enable(gcbf,...
                {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
                'Chord','Chord_units','Xoffset','AeroCenter',...
                'LocationZ','LocationZ_units','modify_smd',...
                'station_delete'});
        end
        % select the new active station
        set(newactive,'Color',[0 0.8 0],'Tag','active station');
        gui_alter(gcbf,'String',num2str(newindex),{'station_current'});
        readActiveStation(newactive);
        if app.modifysmd==true
            cb_modifySMD(cbo,[]);  % first call deselects active station
            cb_modifySMD(cbo,[]);  % second call displays new active
        end
    end
end


function cb_stationnext(cbo,~)
    app = guidata(cbo);  % retrieve the application data
    if app.selectenabled
        % look for active (already selected) station
        active = findobj(gcbf,'Tag','active station');
        if ~isempty(active)
            % station found, remember station number
            index = findCellIndex(app.hg,active);
            % deselect the station
            set(active,'Color',[0 0 0],'Tag','');
            newindex = index + 1;
            if newindex > numel(app.hg)
                newindex = newindex - numel(app.hg);
            end
            newactive = app.hg{newindex};
        else
            newactive = app.hg{1};
            newindex = 1;
            gui_enable(gcbf,...
                {'AirfoilName','TEtype','DegreesTwist','DegreesTwist_units',...
                'Chord','Chord_units','Xoffset','AeroCenter',...
                'LocationZ','LocationZ_units','modify_smd',...
                'station_delete'});
        end
        % select the new active station
        set(newactive,'Color',[0 0.8 0],'Tag','active station');
        gui_alter(gcbf,'String',num2str(newindex),{'station_current'});
        readActiveStation(newactive);
        if app.modifysmd==true
            cb_modifySMD(cbo,[]);  % first call deselects active station
            cb_modifySMD(cbo,[]);  % second call displays new active
        end
    end
end

function cb_modifySMD(cbo,~)
% modify the surface material divisions
    app = guidata(cbo);   % retrieve application data
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    if app.modifysmd==false
        app.modifysmd = true;
        obj = findobj(app.fh,'Tag','modify_smd');
        set(obj,'String','Done modifying');
        gui_disable(gcbf,{'station_new','station_delete','LocationZ','uimenu_rotation'});
        gui_enable(gcbf,{'DPnumber','DPbutton','DPtype','DPsurf','DPchord','DPcdist','DPsdist'});
        gui_enable(gcbf,{'SWnumber','SWbutton','SWmaterial',...
            'SWstation1','SWudp1','SWldp1','SWstation2','SWudp2','SWldp2'});
        set(app.hgGL,'Visible','off');
        cellfun(@(c) set(c,'Visible','off'), app.hg);
        cellfun(@(c) set(c,'Visible','off'), app.swhg);
        
        nDP = numel(app.station(k).dp) - 2;
        uictrl = findobj(app.fh,'Tag','DPnumber');
        dplist = cellstr(num2str(transpose(1:nDP)));
        dplist(nDP+1) = {'add DP'};
        set(uictrl,'String',dplist,'Value',1);
        
        SWindex = [];
        for ksw = 1:numel(app.shearweb)
            if (app.shearweb(ksw).BeginStation==k)
                SWindex(end+1) = ksw;
            end
        end
        nSW = numel(SWindex);
        if nSW > 0
            swlist = cellstr(num2str(transpose(1:nSW)));
        else
            swlist = {''};
        end
        swlist(end+1) = {'add SW'};
        uictrl = findobj(app.fh,'Tag','SWnumber');
        set(uictrl,'String',swlist,'Value',1,'UserData',SWindex);
        
        draw_materialdivisions(app,k);
        readActiveDP(app);
        readActiveSW(app);
    else
        app.modifysmd = false;
        obj = findobj(app.fh,'Tag','modify_smd');
        set(obj,'String','Modify Skin Material Divisions');
        gui_enable(gcbf,{'station_new','station_delete','LocationZ','uimenu_rotation'});
        gui_disable(gcbf,{'DPnumber','DPbutton','DPtype','DPmaterial','DPsurf','DPchord','DPcdist','DPsdist'});
        gui_disable(gcbf,{'SWnumber','SWbutton','SWmaterial',...
            'SWstation1','SWudp1','SWldp1','SWstation2','SWudp2','SWldp2'});
        set(app.hgGL,'Visible','on');
        cellfun(@(c) set(c,'Visible','on'), app.hg);
        cellfun(@(c) set(c,'Visible','on'),  app.swhg);
        cellfun(@(c) set(c,'Visible','off'), app.smdp.hg_sm);
        cellfun(@(c) set(c,'Visible','off'), app.smdp.hg_dp);
        cellfun(@(c) set(c,'Visible','off'), app.smdp.hg_sw);
        legend(app.ax,'off');
        draw_stations(app);
    end
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
end

function readActiveDP(app)
% this function populates the gui controls with the active DP's
% parameter values
    DPnumber = get_uictrl(app.fh,'DPnumber');
    
    % determine which station is active
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    
    % load the stored DP paramater values
    if strcmp(DPnumber.select,'add DP')
        dp.chord = app.newdp.chord;
        dp.type = app.newdp.type;
        dp.material = app.newdp.material;
    else
        dp.chord = app.station(k).dp(1+DPnumber.value);
        dp.type = app.station(k).dptype{1+DPnumber.value};
        dp.material = app.station(k).dpmaterial{1+DPnumber.value};
    end
   
    DPsurf  = findobj(app.fh,'Tag','DPsurf');
    DPchord = findobj(app.fh,'Tag','DPchord');
    DPcdist = findobj(app.fh,'Tag','DPcdist');
    DPsdist = findobj(app.fh,'Tag','DPsdist');
    
    dp.sdist = interp1(app.station(k).c(2:end-1),app.station(k).s(2:end-1),dp.chord)*app.station(k).Chord;
    if dp.chord < 0
        set(DPsurf,'Value',1);  % lower (HP) surface
        dp.chord = abs(dp.chord);
        dp.sdist = -dp.sdist;
    elseif dp.chord > 0
        set(DPsurf,'Value',2);  % upper (LP) surface
    else
        % DP is at leading edge, let user select either HP or LP
    end
    
    DPtype = get_uictrl(app.fh,'DPtype');
    n = strcmp(dp.type,DPtype.string);
    set(DPtype.h,'Value',find(n==1));
    switch dp.type
        case {'flare','hourglass'}
            gui_enable(gcbf,{'DPmaterial'});
        otherwise
            gui_disable(gcbf,{'DPmaterial'});
    end
    
    DPmaterial = get_uictrl(app.fh,'DPmaterial');
    n = strcmp(dp.material,DPmaterial.string);
    if ~any(n)
        n(1) = 1;
    end
    set(DPmaterial.h,'Value',find(n==1));
    
    set(DPchord,'String',num2str(dp.chord*100));  % percent chord
    set(DPcdist,'String',num2str(dp.chord*app.station(k).Chord));  % chordal distance
    set(DPsdist,'String',num2str(dp.sdist));  % surface distance
    
%     app.station(k).x = app.afdb(n).x;
%     app.station(k).y = app.afdb(n).y;
    
end

function readActiveSW(app)
% this function populates the gui controls with the active SW's
% parameter values
    SWnumber = get_uictrl(app.fh,'SWnumber');
    SWmaterial = get_uictrl(app.fh,'SWmaterial');
    SWstation1 = get_uictrl(app.fh,'SWstation1');
    SWstation2 = get_uictrl(app.fh,'SWstation2');
    SWudp1 = get_uictrl(app.fh,'SWudp1');
    SWldp1 = get_uictrl(app.fh,'SWldp1');
    SWudp2 = get_uictrl(app.fh,'SWudp2');
    SWldp2 = get_uictrl(app.fh,'SWldp2');
    
    % determine which station is active
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station

    if k<numel(app.station) && ~isfield(app,'newdp')
        gui_enable(gcbf,{'SWnumber','SWbutton','SWmaterial',...
            'SWstation1','SWudp1','SWldp1','SWstation2','SWudp2','SWldp2'});
    else
        gui_disable(gcbf,{'SWnumber','SWbutton','SWmaterial',...
            'SWstation1','SWudp1','SWldp1','SWstation2','SWudp2','SWldp2'});
        return;
    end
    
    % find upper and lower DPs for the current station and next station
    udp1 = find(app.station( k ).dp(2:end-1)>=0);
    ldp1 = find(app.station( k ).dp(2:end-1)< 0);
    udp2 = find(app.station(k+1).dp(2:end-1)>=0);
    ldp2 = find(app.station(k+1).dp(2:end-1)< 0);

    % update the SW number menu
%     if (SWnumber.value == numel(SWnumber.string)-1) && (SWnumber.value > 1)
%         set(SWnumber.h,'Value',SWnumber.value-1);
%     end
    SWindex = [];
    for ksw = 1:numel(app.shearweb)
        if (app.shearweb(ksw).BeginStation==k)
            SWindex(end+1) = ksw;  % list of shear webs beginning at active station
        end
    end
    nSW = numel(SWindex);
    if nSW > 0
        swlist = cellstr(num2str(transpose(1:nSW)));
    else
        swlist = {''};
    end
    swlist(end+1) = {'add SW'};
%    SWnumber_value = 1;
    if ~isnan(str2double(SWnumber.select))
        SWnumber_value = find(SWindex==SWnumber.userdata(SWnumber.value));
    else
        SWnumber_value = find(strcmp(SWnumber.select,swlist)==1);  % '' or 'add SW'
    end
    set(SWnumber.h,'String',swlist,'Value',SWnumber_value,'UserData',SWindex);
    
    % load the stored SW paramater values
    if strcmp(SWnumber.select,'add SW')
        sw.material = app.newsw.Material;
        sw.begin = app.newsw.BeginStation;
        sw.end = app.newsw.EndStation;
        sw.corner = app.newsw.Corner;
        set(SWmaterial.h,'Value',find(strcmp(sw.material,SWmaterial.string)==1));
        set(SWstation1.h,'String',sw.begin);
        set(SWstation2.h,'String',sw.end);
        set(SWudp1.h,'String',cellstr(num2str(udp1)));
        set(SWldp1.h,'String',cellstr(num2str(ldp1)));
        set(SWudp2.h,'String',cellstr(num2str(udp2)));
        set(SWldp2.h,'String',cellstr(num2str(ldp2)));
        set(SWudp1.h,'Value',find(sw.corner(1)==udp1));
        set(SWldp1.h,'Value',find(sw.corner(2)==ldp1));
        set(SWldp2.h,'Value',find(sw.corner(3)==ldp2));
        set(SWudp2.h,'Value',find(sw.corner(4)==udp2));
    elseif isempty(SWindex)
        set(SWmaterial.h,'Value',1);
        set(SWstation1.h,'String','');
        set(SWstation2.h,'String','');
        set(SWudp1.h,'String',cellstr(num2str(udp1)),'Value',1);
        set(SWldp1.h,'String',cellstr(num2str(ldp1)),'Value',1);
        set(SWudp2.h,'String',cellstr(num2str(udp2)),'Value',1);
        set(SWldp2.h,'String',cellstr(num2str(ldp2)),'Value',1);
    else
        ksw = SWindex(SWnumber_value);
        sw.material = app.shearweb(ksw).Material;
        sw.begin = app.shearweb(ksw).BeginStation;
        sw.end = app.shearweb(ksw).EndStation;
        sw.corner = app.shearweb(ksw).Corner;
        set(SWmaterial.h,'Value',find(strcmp(sw.material,SWmaterial.string)==1));
        set(SWstation1.h,'String',sw.begin);
        set(SWstation2.h,'String',sw.end);
        set(SWudp1.h,'String',cellstr(num2str(udp1)));
        set(SWldp1.h,'String',cellstr(num2str(ldp1)));
        set(SWudp2.h,'String',cellstr(num2str(udp2)));
        set(SWldp2.h,'String',cellstr(num2str(ldp2)));
        set(SWudp1.h,'Value',find(sw.corner(1)==udp1));
        set(SWldp1.h,'Value',find(sw.corner(2)==ldp1));
        set(SWldp2.h,'Value',find(sw.corner(3)==ldp2));
        set(SWudp2.h,'Value',find(sw.corner(4)==udp2));
    end
end

function cb_modifySMD_uictrls(cbo,~)
    app = guidata(cbo);   % retrieve application data
    bview = NuMAD_appdata('get','bview');
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    tag = get(cbo,'Tag');
    switch tag
        case 'DPnumber'
            DPnumber = get_uictrl(app.fh,'DPnumber');
            DPbutton = get_uictrl(app.fh,'DPbutton');
            if strcmp(DPnumber.select,'add DP') 
                if isfield(app,'newdp')
                    return; % do nothing if 'add DP' was already selected
                end
                % modify gui
                set(DPbutton.h,'String','Done');
                gui_alter(gcbf,'BackgroundColor',[1 1 0.8],...
                    {'DPnumber','DPtype','DPmaterial','DPsurf','DPchord',...
                    'DPcdist','DPsdist'})
                gui_disable(gcbf,{'modify_smd','station_previous','station_next','check_model'});
                gui_disable(gcbf,{'SWnumber','SWbutton','SWmaterial',...
                    'SWstation1','SWudp1','SWldp1','SWstation2','SWudp2','SWldp2'});
                % initialize new DP
                app.newdp.chord = 0.01;
                app.newdp.type = 'single';
                app.newdp.material = '';
                app.newdp.hg = line(0,0,0,'Parent',bview.SVhgtf,'Color',[0 0 0],'Marker','o');           
            else
                % the user has selected a DP number
                if isfield(app,'newdp')
                    % stop creating new DP, restore gui to normal
                    set(DPbutton.h,'String','Delete DP');
                    gui_alter(gcbf,'BackgroundColor',[1 1 1],...
                        {'DPnumber','DPtype','DPmaterial','DPsurf','DPchord',...
                        'DPcdist','DPsdist'})
                    gui_enable(gcbf,{'modify_smd','station_previous','station_next','check_model'});
                    gui_enable(gcbf,{'SWnumber','SWbutton','SWmaterial',...
                        'SWstation1','SWudp1','SWldp1','SWstation2','SWudp2','SWldp2'});
                    delete(app.newdp.hg);
                    app = rmfield(app,'newdp');
                end
            end
            
        case 'DPbutton'
            DPnumber = get_uictrl(app.fh,'DPnumber');
            DPbutton = get_uictrl(app.fh,'DPbutton');
            if strcmp(DPnumber.select,'add DP')
                % finish adding a DP - append DP to end
                app.station(k).dp(end+1) = app.newdp.chord;
                app.station(k).dptype{end+1} = app.newdp.type;
                app.station(k).dpmaterial{end+1} = app.newdp.material;
                % sort the DPs
                [app.station(k).dp sortorder] = sort(app.station(k).dp);
                app.station(k).dptype = app.station(k).dptype(sortorder);
                app.station(k).dpmaterial = app.station(k).dpmaterial(sortorder);
                % propagate DP changes into shear web definitons
                for ksw = 1:numel(app.shearweb)
                    if (k==app.shearweb(ksw).BeginStation)
                        corner = 1+app.shearweb(ksw).Corner(1:2);
                        corner(1) = find(corner(1)==sortorder);
                        corner(2) = find(corner(2)==sortorder);
                        app.shearweb(ksw).Corner(1:2) = corner-1;
                    elseif (k==app.shearweb(ksw).EndStation)
                        corner = 1+app.shearweb(ksw).Corner(3:4);
                        corner(1) = find(corner(1)==sortorder);
                        corner(2) = find(corner(2)==sortorder);
                        app.shearweb(ksw).Corner(3:4) = corner-1;
                    end
                end
                % the new DP breaks a segment in two, apply appropriate
                %   skin material to the new piece of that segment
                n = find(sortorder==(DPnumber.value+2))-1;
                sortorder(sortorder==1) = [];
                app.station(k).sm{end+1} = app.station(k).sm{n};
                app.station(k).sm = app.station(k).sm(sortorder-1);
                % update the DP number menu
                set(DPnumber.h,'Value',n);
                nDP = numel(app.station(k).dp) - 2;
                dplist = cellstr(num2str(transpose(1:nDP)));
                dplist(nDP+1) = {'add DP'};
                set(DPnumber.h,'String',dplist);
                delete(app.newdp.hg);
                app = rmfield(app,'newdp');
                % update other gui elements
                set(DPbutton.h,'String','Delete DP');
                gui_alter(gcbf,'BackgroundColor',[1 1 1],...
                    {'DPnumber','DPtype','DPmaterial','DPsurf','DPchord',...
                    'DPcdist','DPsdist'})
                gui_enable(gcbf,{'modify_smd','station_previous','station_next','check_model'});
                gui_enable(gcbf,{'SWnumber','SWbutton','SWmaterial',...
                    'SWstation1','SWudp1','SWldp1','SWstation2','SWudp2','SWldp2'});
            else
                % delete the selected DP
                app.station(k).dp(DPnumber.value+1) = [];
                app.station(k).dptype(DPnumber.value+1) = [];
                app.station(k).dpmaterial(DPnumber.value+1) = [];
                app.station(k).sm(DPnumber.value+1) = [];
                deletelist=[];
                for ksw = 1:numel(app.shearweb)
                    if (k==app.shearweb(ksw).BeginStation)
                        corner = app.shearweb(ksw).Corner(1:2);
                        for kc = 1:2
                            if corner(kc)==(DPnumber.value)
                                deletelist(end+1) = ksw;
                            elseif corner(kc)>(DPnumber.value)
                                app.shearweb(ksw).Corner(kc) = corner(kc)-1;
                            end
                        end
                    elseif (k==app.shearweb(ksw).EndStation)
                        corner = app.shearweb(ksw).Corner(3:4);
                        for kc = 1:2
                            if corner(kc)==(DPnumber.value)
                                deletelist(end+1) = ksw;
                            elseif corner(kc)>(DPnumber.value)
                                app.shearweb(ksw).Corner(kc+2) = corner(kc)-1;
                            end
                        end
                    end
                end
                if ~isempty(deletelist)
                    app.shearweb(deletelist) = [];
                    delete(app.swhg{deletelist});
                    delete(app.smdp.hg_sw{deletelist});
                    app.swhg(deletelist) = [];
                    app.smdp.hg_sw(deletelist) = [];
                end
                
                % update the DP number menu
                if (DPnumber.value == numel(DPnumber.string)-1) && (DPnumber.value > 1)
                    set(DPnumber.h,'Value',DPnumber.value-1)
                end
                nDP = numel(app.station(k).dp) - 2;
                dplist = cellstr(num2str(transpose(1:nDP)));
                dplist(nDP+1) = {'add DP'};
                set(DPnumber.h,'String',dplist);
                if nDP == 0
                    % modify gui
                    set(DPbutton.h,'String','Done');
                    gui_alter(gcbf,'BackgroundColor',[1 1 0.8],...
                        {'DPnumber','DPtype','DPmaterial','DPsurf','DPchord',...
                        'DPcdist','DPsdist'})
                    % initialize new DP
                    app.newdp.chord = 0;
                    app.newdp.type = 'single';
                    app.newdp.hg = line(0,0,0,'Parent',app.ax,'Color',[0 0 0],'Marker','o');
                end
            end
            
        case 'DPtype'
            DPnumber = get_uictrl(app.fh,'DPnumber');
            DPtype = get_uictrl(app.fh,'DPtype');
            switch DPtype.select
                case {'flare','hourglass'}
                    if strcmp(DPnumber.select,'add DP')
                        app.newdp.type = DPtype.select;
                        app.newdp.material = '**UNSPECIFIED**';
                    else
                        app.station(k).dptype{DPnumber.value+1} = DPtype.select;
                        app.station(k).dpmaterial{DPnumber.value+1} = '**UNSPECIFIED**';
                    end
                    gui_enable(gcbf,{'DPmaterial'});
                otherwise
                    if strcmp(DPnumber.select,'add DP')
                        app.newdp.type = DPtype.select;
                        app.newdp.material = '';
                    else
                        app.station(k).dptype{DPnumber.value+1} = DPtype.select;
                        app.station(k).dpmaterial{DPnumber.value+1} = '';
                    end
                    gui_disable(gcbf,{'DPmaterial'});
            end
            
        case 'DPmaterial'
            DPnumber = get_uictrl(app.fh,'DPnumber');
            DPmaterial = get_uictrl(app.fh,'DPmaterial');
            if strcmp(DPnumber.select,'add DP')
                app.newdp.material = DPmaterial.select;
            else
                app.station(k).dpmaterial{DPnumber.value+1} = DPmaterial.select;
            end
            
        case 'DPsurf'
            DPnumber = get_uictrl(app.fh,'DPnumber');
            DPchord = get_uictrl(app.fh,'DPchord');
            DPsurf = get_uictrl(app.fh,'DPsurf');
            
            if strcmp(DPnumber.select,'add DP')
                dp.chord = app.newdp.chord;  % get current stored value of new DP
            else
                dp.chord = app.station(k).dp(DPnumber.value+1);  % get current stored value
            end
            
            if ((dp.chord < 0) && (DPsurf.value==2)) || ((dp.chord > 0) && (DPsurf.value==1))
                dp.chord = -1 * dp.chord;  % switch surface
            end
            
            if strcmp(DPnumber.select,'add DP')
                app.newdp.chord = dp.chord;
            else
                % update the chord value
                app.station(k).dp(DPnumber.value+1) = dp.chord;
                % re-sort the DPs
                [app.station(k).dp sortorder] = sort(app.station(k).dp);
                app.station(k).dptype = app.station(k).dptype(sortorder);
                set(DPnumber.h,'Value',find(sortorder==(DPnumber.value+1))-1);
            end
            
        case 'DPchord'
            DPnumber = get_uictrl(app.fh,'DPnumber');
            DPchord = get_uictrl(app.fh,'DPchord');
            DPsurf = get_uictrl(app.fh,'DPsurf');
            dp.chord = app.station(k).dp(DPnumber.value+1);  % get current stored value
            
            num = str2double(DPchord.string);  % try to convert to number
            if isnan(num)
                % input was not valid => revert to stored value
                str = sprintf('%g',dp.chord);
                set(DPchord.h,'string',str);
            else
                % input valid => save new value
                if DPsurf.value==2  % upper (LP) surface
                    dp.chord = abs(num)/100;
                else  % lower (HP) surface
                    dp.chord = -abs(num)/100;
                end
                if strcmp(DPnumber.select,'add DP')
                    app.newdp.chord = dp.chord;
                else
                    % update the chord value
                    app.station(k).dp(DPnumber.value+1) = dp.chord;
                    % re-sort the DPs
                    [app.station(k).dp sortorder] = sort(app.station(k).dp);
                    app.station(k).dptype = app.station(k).dptype(sortorder);
                    set(DPnumber.h,'Value',find(sortorder==(DPnumber.value+1))-1);
                end
            end
            
        case 'DPcdist'
            DPnumber = get_uictrl(app.fh,'DPnumber');
            DPcdist = get_uictrl(app.fh,'DPcdist');
            DPsurf = get_uictrl(app.fh,'DPsurf');
            dp.chord = app.station(k).dp(DPnumber.value+1);  % get current stored value
            
            num = str2double(DPcdist.string);  % try to convert to number
            if isnan(num)
                % input was not valid => revert to stored value
                str = sprintf('%g',dp.chord*app.station(k).Chord);
                set(DPcdist.h,'string',str);
            else
                % input valid => save new value
                if DPsurf.value==2  % upper (LP) surface
                    dp.chord = abs(num)/app.station(k).Chord;
                else  % lower (HP) surface
                    dp.chord = -abs(num)/app.station(k).Chord;
                end
                if strcmp(DPnumber.select,'add DP')
                    app.newdp.chord = dp.chord;
                else
                    app.station(k).dp(DPnumber.value+1) = dp.chord;
                end
            end
            
        case 'DPsdist'
            DPnumber = get_uictrl(app.fh,'DPnumber');
            DPsdist = get_uictrl(app.fh,'DPsdist');
            DPsurf = get_uictrl(app.fh,'DPsurf');
            dp.chord = app.station(k).dp(DPnumber.value+1);  % get current stored value
            
            num = str2double(DPsdist.string);  % try to convert to number
            if isnan(num)
                % input was not valid => revert to stored value
                dp.sdist = interp1(app.station(k).c(2:end-1),app.station(k).s(2:end-1),dp.chord)*app.station(k).Chord;
                if dp.chord < 0
                    dp.sdist = -dp.sdist;
                end
                str = sprintf('%g',abs(dp.sdist));
                set(DPsdist.h,'string',str);
            else
                % input valid => save new value
                if DPsurf.value==2  % upper (LP) surface
                    dp.chord = interp1(app.station(k).s(2:end-1),app.station(k).c(2:end-1), abs(num)/app.station(k).Chord);
                else  % lower (HP) surface
                    dp.chord = interp1(app.station(k).s(2:end-1),app.station(k).c(2:end-1),-abs(num)/app.station(k).Chord);
                end
                if strcmp(DPnumber.select,'add DP')
                    app.newdp.chord = dp.chord;
                else
                    app.station(k).dp(DPnumber.value+1) = dp.chord;
                end
            end
            
    end
    
    % refresh the gui
    readActiveDP(app);
    readActiveSW(app);  % also need to refresh SW because DP numbers can change
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
    
    draw_materialdivisions(app,k)
end

function cb_modifySW_uictrls(cbo,~)
    app = guidata(cbo);   % retrieve application data
    bview = NuMAD_appdata('get','bview');
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    
    SWnumber = get_uictrl(app.fh,'SWnumber');
    SWbutton = get_uictrl(app.fh,'SWbutton');
    SWmaterial = get_uictrl(app.fh,'SWmaterial');
    SWstation1 = get_uictrl(app.fh,'SWstation1');
    SWstation2 = get_uictrl(app.fh,'SWstation2');
    SWudp1 = get_uictrl(app.fh,'SWudp1');
    SWldp1 = get_uictrl(app.fh,'SWldp1');
    SWudp2 = get_uictrl(app.fh,'SWudp2');
    SWldp2 = get_uictrl(app.fh,'SWldp2');
    
    tag = get(cbo,'Tag');
    switch tag
        case 'SWnumber'
            if strcmp(SWnumber.select,'add SW')
                if isfield(app,'newsw')
                    return; % do nothing if 'add SW' was already selected
                end
                % modify gui
                set(SWbutton.h,'String','Done');
                gui_alter(gcbf,'BackgroundColor',[1 1 0.8],...
                    {'SWnumber','SWmaterial','SWudp1','SWldp1',...
                     'SWudp2','SWldp2'})
                gui_disable(gcbf,{'modify_smd','station_previous','station_next','check_model'});
                gui_disable(gcbf,{'DPnumber','DPbutton','DPtype','DPmaterial','DPsurf','DPchord','DPcdist','DPsdist'}); 
                % initialize new SW
                app.newsw.Material = '**UNSPECIFIED**';
                app.newsw.BeginStation = k;
                app.newsw.EndStation = k+1;
                app.newsw.Corner = str2double([SWudp1.string(1) SWldp1.string(end) SWldp2.string(end) SWudp2.string(1)]);
                app.hg_newsw = line(0,0,0,'Parent',bview.SVhgtf,'Color',[1 0 0],'Visible','off');
            else
                % the user has selected a SW number
                if isfield(app,'newsw')
                    % stop creating new SW, restore gui to normal
                    set(SWbutton.h,'String','Delete SW');
                    gui_alter(gcbf,'BackgroundColor',[1 1 1],...
                        {'SWnumber','SWmaterial','SWudp1','SWldp1',...
                         'SWudp2','SWldp2'})
                    gui_enable(gcbf,{'modify_smd','station_previous','station_next','check_model'});
                    gui_enable(gcbf,{'DPnumber','DPbutton','DPtype','DPmaterial','DPsurf','DPchord','DPcdist','DPsdist'});
                    delete(app.hg_newsw);
                    app = rmfield(app,'newsw');
                end
            end
            
        case 'SWbutton'
            if strcmp(SWnumber.select,'add SW')
                % the user has finished creating a new SW
                set(SWbutton.h,'String','Delete SW');
                gui_alter(gcbf,'BackgroundColor',[1 1 1],...
                        {'SWnumber','SWmaterial','SWudp1','SWldp1',...
                         'SWudp2','SWldp2'})
                gui_enable(gcbf,{'modify_smd','station_previous','station_next','check_model'});
                gui_enable(gcbf,{'DPnumber','DPbutton','DPtype','DPmaterial','DPsurf','DPchord','DPcdist','DPsdist'});
                % append SW to end
                if isempty(app.shearweb)
                    app.shearweb = app.newsw;
                else
                    app.shearweb(end+1) = app.newsw;
                end
                app.swhg{end+1} = line(0,0,0,'Parent',bview.hgtf,...
                    'Color',[0 0 1],'LineWidth',1.0);
                app.smdp.hg_sw{end+1} = line(0,0,0,'Parent',bview.SVhgtf,...
                    'Color',[0 0 1],'LineWidth',1.0);
                % clean up
                delete(app.hg_newsw);
                app = rmfield(app,'newsw');
                
                % update SW number menu
                SWindex = [];
                for ksw = 1:numel(app.shearweb)
                    if (app.shearweb(ksw).BeginStation==k)
                        SWindex(end+1) = ksw;
                    end
                end
                nSW = numel(SWindex);
                if nSW > 0
                    swlist = cellstr(num2str(transpose(1:nSW)));
                else
                    swlist = {''};
                end
                swlist(end+1) = {'add SW'};
                uictrl = findobj(app.fh,'Tag','SWnumber');
                set(uictrl,'String',swlist,'Value',1,'UserData',SWindex);
                
            else
                % the user wants to delete a SW
                ksw = SWnumber.userdata(SWnumber.value);
                delete(app.swhg{ksw});
                delete(app.smdp.hg_sw{ksw});
                app.swhg(ksw) = [];
                app.smdp.hg_sw(ksw) = [];
                app.shearweb(ksw) = [];
                
                % update SW number menu
                SWindex = [];
                for ksw = 1:numel(app.shearweb)
                    if (app.shearweb(ksw).BeginStation==k)
                        SWindex(end+1) = ksw;
                    end
                end
                nSW = numel(SWindex);
                if nSW > 0
                    swlist = cellstr(num2str(transpose(1:nSW)));
                else
                    swlist = {''};
                end
                swlist(end+1) = {'add SW'};
                uictrl = findobj(app.fh,'Tag','SWnumber');
                set(uictrl,'String',swlist,'Value',1,'UserData',SWindex);
            end
            
        case 'SWmaterial'
            if strcmp(SWnumber.select,'add SW')
                app.newsw.Material = SWmaterial.select;
            else
                ksw = SWnumber.userdata(SWnumber.value);  % active shear web
                app.shearweb(ksw).Material = SWmaterial.select;
            end
            
        case {'SWudp1','SWldp1','SWldp2','SWudp2'}
            n = strcmp(tag,{'SWudp1','SWldp1','SWldp2','SWudp2'});
            uictrl = get_uictrl(app.fh,tag);
            if strcmp(SWnumber.select,'add SW')
                app.newsw.Corner(n) = str2double(uictrl.select);
            else
                ksw = SWnumber.userdata(SWnumber.value);  % active shear web
                app.shearweb(ksw).Corner(n) = str2double(uictrl.select);
            end
    end
    
    % refresh the gui
    readActiveSW(app)
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
    
    draw_materialdivisions(app,k)
end    


function cb_smmodify(cbo,~)
    %jcb: The cbo passed to this function corresponds to the context menu
    %     item and not the graphics object.  Need to save the current
    %     object right away before user has time to click on something else
    obj = gco;  % save the handle of the current object
    
    app = guidata(cbo);   % retrieve application data
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);  % get index of active station
    n = findCellIndex(app.smdp.hg_sm(2,:),obj);  % get index of surface segment
    newmaterial = get(cbo,'Label');  % get name of material chosen from context menu label
    
    %c = strcmp(newmaterial,app.active.list);  % get the index of the material color
    %set(obj,'Color',app.active.color{c});
    app.station(k).sm{n} = newmaterial;
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
    
    draw_materialdivisions(app,k)
end

function cb_dpmodify(cbo,~)
    %jcb: The cbo passed to this function corresponds to the context menu
    %     item and not the graphics object.  Need to save the current
    %     object right away before user has time to click on something else
    obj = gco;  % save the handle of the current object
    app = guidata(cbo);   % retrieve application data
    active = findobj(gcbf,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);    % get index of active station
    n = findCellIndex(app.smdp.hg_dp(2,:),obj);  % get index of delineation point
    label = get(cbo,'Label');  % get the menu item chosen
    switch label
        case 'Single'
            app.station(k).dptype{n} = 'single';
            set(obj,'Marker','o','MarkerSize',6);
        case 'Double'
            app.station(k).dptype{n} = 'double';
            set(obj,'Marker','v','MarkerSize',7);
        case 'Hourglass'
            app.station(k).dptype{n} = 'hourglass';
            set(obj,'Marker','d','MarkerSize',7);
        case 'Flare'
            app.station(k).dptype{n} = 'flare';
            set(obj,'Marker','^','MarkerSize',7);
    end
    DPnumber = findobj(app.fh,'Tag','DPnumber');
    set(DPnumber,'Value',n-1);
    readActiveDP(app);
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
end

function cb_dpmove(cbo,~)
    %jcb: The cbo passed to this function corresponds to the context menu
    %     item and not the graphics object.  Need to save the current
    %     object right away before user has time to click on something else
    obj = gco;  % save the handle of the current object
    app = guidata(cbo);   % retrieve application data
    active = findobj(app.fh,'Tag','active station');  % get handle of active station
    k = findCellIndex(app.hg,active);    % get index of active station
    n = findCellIndex(app.smdp.hg_dp(2,:),obj);  % get index of delineation point
    label = get(cbo,'Label');  % get the menu item chosen

    sta = app.station(k);
    x = (sta.xy(:,1) - sta.Xoffset) * sta.Chord;
    y = (sta.xy(:,2)              ) * sta.Chord;
    twist = sta.DegreesTwist * pi/180;
    xt = cos(twist) * x - sin(twist) * y;
    yt = sin(twist) * x + cos(twist) * y;
    %zt = ones(size(x)) * sta.LocationZ;

%     % find the normalized position of each point
%     xn = sta.x;
%     [~,j] = min(xn);
%     if sta.y(j) >= 0
%         j = j-1;
%     end
%     xn(1:j) = -1*xn(1:j);   
    xn = sta.c;
    
    app.motion.staID = k;
    app.motion.dpID = n;
    app.motion.dp = sta.dp(n);
    app.motion.xn = xn;
    app.motion.xt = xt;
    app.motion.yt = yt;
    app.motion.hg = line(0,0,0,'Parent',app.ax,'Color',[0 0 0],'Marker','o');
    
    pointer = get(gca,'CurrentPoint');
    app.motion.datum = pointer(1,:);
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
    set(app.fh,'WindowButtonMotionFcn',@wbmf_dpmove);
    set(app.fh,'WindowButtonDownFcn',@wbdf_dpmovefinish);

end

function wbmf_dpmove(cbo,~)
% WindowButtonMotionFcn
    app = guidata(cbo);   % retrieve application data
    % get current position of mouse pointer
    pointer = get(gca,'CurrentPoint');  
    % calculate a number that represents the amount of pointer movement
    movement = sum((app.motion.datum - pointer(1,:))); 
    % use this number to define dp location
    newdp = app.motion.dp + round(movement*100)/100;  
    if newdp < -0.99, newdp = -0.99; end
    if newdp >  0.99, newdp =  0.99; end
    xt = interp1(app.motion.xn,app.motion.xt,newdp);
    yt = interp1(app.motion.xn,app.motion.yt,newdp);
    set(app.motion.hg,'YData',xt,...
                      'ZData',yt);
    %fprintf('%f\n',movement);
    
%     DPnumber = findobj(app.fh,'Tag','DPnumber');
%     set(DPnumber,'Value',app.motion.dpID-1);
%     readActiveDP(app);
end

function wbdf_dpmovefinish(cbo,~)
% WindowButtonDownFcn
    app = guidata(cbo);   % retrieve application data
    % get current position of mouse pointer
    pointer = get(gca,'CurrentPoint');  
    % calculate a number that represents the amount of pointer movement
    movement = sum((app.motion.datum - pointer(1,:))); 
    % use this number to define dp location
    newdp = app.motion.dp + round(movement*100)/100;  
    if newdp < -0.99, newdp = -0.99; end
    if newdp >  0.99, newdp =  0.99; end
    
    % update the stored information
    app.station(app.motion.staID).dp(app.motion.dpID) = newdp;
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end

    DPnumber = findobj(app.fh,'Tag','DPnumber');
    set(DPnumber,'Value',app.motion.dpID-1);
    readActiveDP(app);
    
    % cleanup and redraw
    set(app.fh,'WindowButtonMotionFcn','');
    set(app.fh,'WindowButtonDownFcn','');
    delete(app.motion.hg);
    draw_materialdivisions(app,app.motion.staID);
end


function angle = radwrap(angle)
% wrap radian angles to -pi..+pi limits
    if angle > pi
        angle = angle - 2*pi;
    elseif angle < -pi
        angle = angle + 2*pi;
    end
end


function wmf_pan(hObject,eventdata)
    currpt = get(gca,'CurrentPoint'); % get mouse point
    bview = NuMAD_appdata('get','bview');
    
    Raz = makehgtform('zrotate',bview.az);
    Rel = makehgtform('xrotate',bview.el);
    R = Rel*Raz;
    
   warnstate = warning('off','MATLAB:nearlySingularMatrix');
   warning('off','MATLAB:singularMatrix');
    mp = R([1,2],[1,2])\currpt(1,1:2)';
    mpf = R([1,2],[1,2])\bview.firstpt(1,1:2)';
    xyz{1} = [mpf(1) mp(1); mpf(2) mp(2); 0 0; 1 1];
    xyznorm(1) = 0.5*sum((xyz{1}(:,1)-xyz{1}(:,2)).^2);
    
    mp = R([1,2],[1,3])\currpt(1,1:2)';
    mpf = R([1,2],[1,3])\bview.firstpt(1,1:2)';
    xyz{2} = [mpf(1) mp(1); 0 0; mpf(2) mp(2); 1 1];
    xyznorm(2) = sum((xyz{2}(:,1)-xyz{2}(:,2)).^2);
    
    mp = R([1,2],[2,3])\currpt(1,1:2)';
    mpf = R([1,2],[2,3])\bview.firstpt(1,1:2)';
    xyz{3} = [0 0; mpf(1) mp(1); mpf(2) mp(2); 1 1];
    xyznorm(3) = sum((xyz{2}(:,1)-xyz{2}(:,2)).^2);
   warning(warnstate);
    
    [~,ki] = min(xyznorm);
    bview.movement = xyz{ki}(1:3,2)-xyz{ki}(1:3,1);
    if ki~=1 % distinguish planar and vertical motions
        bview.movement(1:2) = [0;0];
    end
    
    if isequal(bview.modifysmd,false)
        trans = bview.trans + bview.movement;
        SVtrans = bview.SVtrans;
    else
        trans = bview.trans;
        SVtrans = bview.SVtrans + bview.movement;
    end
    T = makehgtform('translate',trans(1),trans(2),trans(3));
    S = makehgtform('scale',bview.scale);
    set(bview.hgtf,'Matrix',Rel*Raz*T*S);
    T = makehgtform('translate',SVtrans(1),SVtrans(2),SVtrans(3));
    S = makehgtform('scale',bview.SVscale);
    set(bview.SVhgtf,'Matrix',Rel*Raz*T*S);
    
    NuMAD_appdata('set','bview',bview);
end

function wmf_orbit(hObject,eventdata)
    P1 = get(gcbf,'CurrentPoint');
    bview = NuMAD_appdata('get','bview');
    motionvector = P1-bview.P0;

    bview.az = radwrap(bview.az + 0.005*motionvector(1));
    bview.el = radwrap(bview.el - 0.002*motionvector(2));

    Raz = makehgtform('zrotate',bview.az);
    Rel = makehgtform('xrotate',bview.el);
    T = makehgtform('translate',bview.trans(1),bview.trans(2),bview.trans(3));
    S = makehgtform('scale',bview.scale);
    set(bview.hgtf,'Matrix',Rel*Raz*T*S);
    T = makehgtform('translate',bview.SVtrans(1),bview.SVtrans(2),bview.SVtrans(3));
    S = makehgtform('scale',bview.SVscale);
    set(bview.SVhgtf,'Matrix',Rel*Raz*T*S);
    bview.P0=P1;
    NuMAD_appdata('set','bview',bview);
end

function wmf_zoom(hObject,eventdata)
    P1 = get(gcbf,'CurrentPoint');
    bview = NuMAD_appdata('get','bview');
    motionvector = P1-bview.P0;
    bview.P0=P1;
    
    multiplier = (1-0.01*motionvector(2));
    
    Raz = makehgtform('zrotate',bview.az);
    Rel = makehgtform('xrotate',bview.el);
    if isequal(bview.modifysmd,false)
        bview.scale = bview.scale*multiplier;
        bview.trans = bview.trans*multiplier;
        if bview.scale <= 0
            bview.scale = eps;
        end
        T = makehgtform('translate',bview.trans(1),bview.trans(2),bview.trans(3));
        S = makehgtform('scale',bview.scale);
        set(bview.hgtf,'Matrix',Rel*Raz*T*S);
    else
        bview.SVscale = bview.SVscale*multiplier;
        bview.SVtrans = bview.SVtrans*multiplier;
        if bview.SVscale <= 0
            bview.SVscale = eps;
        end
        T = makehgtform('translate',bview.SVtrans(1),bview.SVtrans(2),bview.SVtrans(3));
        S = makehgtform('scale',bview.SVscale);
        set(bview.SVhgtf,'Matrix',Rel*Raz*T*S);
    end
    NuMAD_appdata('set','bview',bview);
end
    

function wf_scrollwheel(hObject,eventdata)
    bview = NuMAD_appdata('get','bview');
    Raz = makehgtform('zrotate',bview.az);
    Rel = makehgtform('xrotate',bview.el);
    if isequal(bview.modifysmd,false)
        bview.scale = bview.scale*(1 + 0.1*eventdata.VerticalScrollCount);
        bview.trans = bview.trans*(1 + 0.1*eventdata.VerticalScrollCount);
        if bview.scale <= 0
            bview.scale = eps;
        end
        T = makehgtform('translate',bview.trans(1),bview.trans(2),bview.trans(3));
        S = makehgtform('scale',bview.scale);
        set(bview.hgtf,'Matrix',Rel*Raz*T*S);
    else
        bview.SVscale = bview.SVscale*(1 + 0.1*eventdata.VerticalScrollCount);
        bview.SVtrans = bview.SVtrans*(1 + 0.1*eventdata.VerticalScrollCount);
        if bview.SVscale <= 0
            bview.SVscale = eps;
        end
        T = makehgtform('translate',bview.SVtrans(1),bview.SVtrans(2),bview.SVtrans(3));
        S = makehgtform('scale',bview.SVscale);
        set(bview.SVhgtf,'Matrix',Rel*Raz*T*S);
    end
    NuMAD_appdata('set','bview',bview);
end

function wf_buttondown(hObject,eventdata)
    app = guidata(hObject);   % retrieve application data
    bview = NuMAD_appdata('get','bview');
    bview.modifysmd = app.modifysmd;
    selType = get(gcbf,'SelectionType');
    motionfcn = '';
    switch lower(selType);
        case 'normal'
            bview.lastclick = 'normal';   
            motionfcn = 'translate';
        case 'extend'
            bview.lastclick = 'extend';
            motionfcn = 'rotate';
        case 'open'
            if isequal(bview.lastclick,'extend')
                bview.lastclick = '';
                motionfcn = 'reset';
            end
        otherwise
            bview.lastclick = '';
    end
    if ~isempty(bview.function) 
        if isequal(lower(selType),'normal')
            motionfcn = bview.function;
        end
    end
    switch motionfcn
        case 'translate'
            bview.firstpt = get(gca,'CurrentPoint'); % get mouse point 
            bview.movement = zeros(3,1);
            set(gcbf, 'WindowButtonMotionFcn', @wmf_pan); 
        case 'rotate'
            bview.P0 = get(gcbf,'CurrentPoint');
            set(gcbf, 'WindowButtonMotionFcn', @wmf_orbit);
        case 'zoom'
            bview.P0 = get(gcbf,'CurrentPoint');
            set(gcbf, 'WindowButtonMotionFcn', @wmf_zoom);
        case 'reset'
                bview.az = -0.8;
                bview.el = -0.9;
                bview.scale = 250/(1 + bview.bladetip);
                bview.SVscale = bview.scale;
                midLocZ = 0.5*bview.bladetip;
                bview.trans = [-midLocZ*bview.scale; 0; 0];
                bview.SVtrans = zeros(3,1);
                bview.movement = zeros(3,1);
                Raz = makehgtform('zrotate',bview.az);
                Rel = makehgtform('xrotate',bview.el);
                T = makehgtform('translate',bview.trans(1),bview.trans(2),bview.trans(3));
                S = makehgtform('scale',bview.scale);
                set(bview.hgtf,'Matrix',Rel*Raz*T*S);
                T = makehgtform('translate',bview.SVtrans(1),bview.SVtrans(2),bview.SVtrans(3));
                S = makehgtform('scale',bview.SVscale);
                set(bview.SVhgtf,'Matrix',Rel*Raz*T*S);
    end
    NuMAD_appdata('set','bview',bview);
end

function wf_buttonup(hObject,eventdata)
    motionfcn = get(gcbf, 'WindowButtonMotionFcn');
    set(gcbf, 'WindowButtonMotionFcn', '' );
    bview = NuMAD_appdata('get','bview');
    if isequal(motionfcn,@wmf_pan)
        if isequal(bview.modifysmd,false)
            bview.trans = bview.trans + bview.movement;
        else
            bview.SVtrans = bview.SVtrans + bview.movement;
        end
        NuMAD_appdata('set','bview',bview);
    end
end

function wf_keypress(hObject,eventdata)

end

function wf_keyrelease(hObject,eventdata)
    if any(strcmp(eventdata.Key,{'r','t','g','z'}))
        bview = NuMAD_appdata('get','bview');
        switch eventdata.Key
            case {'t','g'}
                if isequal(bview.function,'translate')
                    bview.function = '';  %toggle off
                else
                    bview.function = ''; %toggle on (translate default)
                end
            case 'r'
                if isequal(bview.function,'rotate')
                    bview.function = ''; %toggle off
                else
                    bview.function = 'rotate'; %toggle rotate on
                end
            case 'z'
                if isequal(bview.function,'zoom')
                    bview.function = ''; %toggle off
                else
                    bview.function = 'zoom'; %toggle zoom on
                end
        end
        hLabels = findobj(gcbf,'tag','uimenu_viewtool');
        set(hLabels,'checked','off');
        switch bview.function
            case ''
                set(gcbf,'Pointer','arrow');
                set(hLabels(3),'checked','on');
            case 'rotate'
                set(gcf,'Pointer','custom','PointerShapeCData',bview.pointer.rotate);
                set(hLabels(2),'checked','on');
            case 'zoom'
                set(gcf,'Pointer','custom','PointerShapeCData',bview.pointer.zoom);
                set(hLabels(1),'checked','on');
        end     
        NuMAD_appdata('set','bview',bview);
    end
end

function cb_viewtool(hObject,eventdata)
    label = get(hObject,'Label');
    hLabels = findobj(gcbf,'tag','uimenu_viewtool');
    set(hLabels,'checked','off');
    set(hObject,'checked','on');
    bview = NuMAD_appdata('get','bview');
    switch label
        case 'Translate tool [T]'
            bview.function = ''; %translate is default view tool function
            set(gcbf,'Pointer','arrow');
        case 'Rotate 3D tool [R]'
            bview.function = 'rotate'; %toggle rotate on
            set(gcf,'Pointer','custom','PointerShapeCData',bview.pointer.rotate);
        case 'Zoom tool [Z]'
            bview.function = 'zoom'; %toggle zoom on
            set(gcf,'Pointer','custom','PointerShapeCData',bview.pointer.zoom);
    end
    NuMAD_appdata('set','bview',bview);
end

function cb_orientation(hObject,eventdata)
    blade = NuMAD_appdata('get','blade');
    bview = NuMAD_appdata('get','bview');
    label = get(hObject,'Label');
    switch label
        case 'View Upper (LP)'
            bview.az = 0;
            bview.el = 0;
        case 'View Lower (HP)'
            bview.az = 0;
            bview.el = pi;
        case 'View L.E.'
            if isequal(blade.BladeRotation,'cw')
                bview.az = pi;
                bview.el = -pi/2;
            else
                bview.az = 0;
                bview.el = -pi/2;
            end
        case 'View T.E.'
            if isequal(blade.BladeRotation,'cw')
                bview.az = 0;
                bview.el = -pi/2;
            else
                bview.az = pi;
                bview.el = -pi/2;
            end
        case 'View Root'
            bview.az = pi/2;
            bview.el = -pi/2;
        case 'View Tip'
            bview.az = -pi/2;
            bview.el = -pi/2;
        case 'Reset View'
            bview.az = -0.8;
            bview.el = -0.9;
            bview.scale = 250/(1 + bview.bladetip);
            bview.SVscale = bview.scale;
            midLocZ = 0.5*bview.bladetip;
            bview.trans = [-midLocZ*bview.scale; 0; 0];
            bview.SVtrans = zeros(3,1);
            bview.movement = zeros(3,1);
    end
    Raz = makehgtform('zrotate',bview.az);
    Rel = makehgtform('xrotate',bview.el);
    T = makehgtform('translate',bview.trans(1),bview.trans(2),bview.trans(3));
    S = makehgtform('scale',bview.scale);
    set(bview.hgtf,'Matrix',Rel*Raz*T*S);
    T = makehgtform('translate',bview.SVtrans(1),bview.SVtrans(2),bview.SVtrans(3));
    S = makehgtform('scale',bview.SVscale);
    set(bview.SVhgtf,'Matrix',Rel*Raz*T*S);
    NuMAD_appdata('set','bview',bview);
end

function update_activelist(fh,active)
    app = guidata(fh);  % retrieve application data
    app.active = active;
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
    
    delete(get(app.sm_cmenu,'Children'));
    for k = 1:numel(app.active.list)
        uimenu(app.sm_cmenu,'Label',app.active.list{k},...
                            'Callback',@cb_smmodify);
    end
    
%     % populate the list of active materials
%     uictrl = get_uictrl(app.fh,'DPmaterial');
%     activelist = unique([app.active.list; uictrl.select]);
%     value = find(strcmp(app.active.list,uictrl.select)==1);
%     set(uictrl.h,'String',activelist,'Value',value);
%     
%     uictrl = get_uictrl(app.fh,'SWmaterial');
%     activelist = unique([app.active.list; uictrl.select]);
%     value = find(strcmp(app.active.list,uictrl.select)==1);
%     set(uictrl.h,'String',activelist,'Value',value);
end

function update_paths(fh,pathconfig)
    app = guidata(fh);  % retrieve application data
    app.settings.ansys_path = pathconfig.ansys_path;
    app.settings.ansys_product = pathconfig.ansys_product;
    app.settings.bmodes_path = pathconfig.bmodes_path;
    app.settings.precomp_path = pathconfig.precomp_path;
    settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
    settings.ansys_path = pathconfig.ansys_path;
    settings.ansys_product = pathconfig.ansys_product;
    settings.bmodes_path = pathconfig.bmodes_path;
    settings.precomp_path = pathconfig.precomp_path;
    writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
end

function update_ansys(fh,ansys)
    app = guidata(fh);  % retrieve application data
    app.ansys = ansys;
    app.settings.ansys_path = ansys.path;
    app.settings.ansys_product = ansys.product;
    settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
    settings.ansys_path = ansys.path;
    settings.ansys_product = ansys.product;
    writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
end

function cb_close(cbo,~)
    app = guidata(cbo);  % load application data
    try
        position = get(app.fh,'Position');
        settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
        settings.xy_main = position;
        if isdir(app.settings.job_path)
            settings.job_path = app.settings.job_path;
        end
        % CK prevent parallel write error
        t = getCurrentTask;
        if isempty(t)
            writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
        end
        NuMAD_appdata('cleanup');
    catch ME
        h=warndlg('Could not complete close request.');
        uiwait(h);
        try
            NuMAD_appdata('cleanup');
        catch ME2
            disp('NuMAD_appdata cleanup error:');
            disp(ME2)
        end
        closereq;
        rethrow(ME);
    end
    closereq;
end

function cb_bladerotation(cbo,~)
    app = guidata(cbo);  % load application data
    blade = NuMAD_appdata('get','blade');
    
    switch get(cbo,'tag')
        case 'uimenu_cw'
            app.BladeRotation = 'cw';
            blade.BladeRotation = 'cw';
            gui_alter(app.fh,'checked','on',{'uimenu_cw'});
            gui_alter(app.fh,'checked','off',{'uimenu_ccw'});
            %view(app.ax,135,30);
        case 'uimenu_ccw'
            app.BladeRotation = 'ccw';
            blade.BladeRotation = 'ccw';
            gui_alter(app.fh,'checked','off',{'uimenu_cw'});
            gui_alter(app.fh,'checked','on',{'uimenu_ccw'});
            %view(app.ax,45,30);
    end
    
    NuMAD_appdata('set','blade',blade);
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
    
    draw_stations(app);
end

function cb_generate(cbo,~)
    app = guidata(cbo);  % load application data
    blade = NuMAD_appdata('get','blade');
    
    if isempty(app.settings.job_name)
        warndlg('Please save the model before generating output.','Operation Not Permitted');
        return;
    end
    cb_checkbladedata(cbo,[]);
    cb_checkmatdata(cbo,[]);
    app = guidata(cbo);  % load application data again after cb_checkbladedata and cb_checkmatdata
    
    if app.checkpassed
        if app.ansys.shell7gen
            shell7_name = 'shell7.src';
        else
            shell7_name = '_TMP_shell7.src';
        end
        filename = fullfile(app.settings.job_path,shell7_name);
        %write_shell7(app,blade,filename);
        if app.ansys.dbgen
            if isempty(app.settings.ansys_path)
                errordlg('Path to ANSYS not specified. Aborting.','Operation Not Permitted');
                return;
            end
            old_dir = pwd;
            try
                %tcl: exec "$ANSYS_path" -b -p $AnsysProductVariable -I shell7.src -o output.txt
                cd(app.settings.job_path);
                ansys_call = sprintf('"%s" -b -p %s -I %s -o output.txt',...
                    app.settings.ansys_path,app.settings.ansys_product,shell7_name);
                %[status,result] = dos(ansys_call);  % the windows system call to run the above ansys command
                
%                 if isequal(status,0)
%                     % dos command completed successfully; log written to output.txt
%                     if app.batchrun
%                         disp('ANSYS batch run to generate database (.db) has completed. See "output.txt" for any warnings.');
%                     else
%                         helpdlg('ANSYS batch run to generate database (.db) has completed. See "output.txt" for any warnings.','ANSYS Call Completed');
%                     end
%                 end
%                 if isequal(status,7)
%                     % an error has occured which is stored in output.txt 
%                     if app.batchrun
%                         disp('Could not complete ANSYS call. See "output.txt" for details.');
%                     else
%                         warndlg('Could not complete ANSYS call. See "output.txt" for details.','Error: ANSYS Call');
%                     end
%                 end
                if 0==app.ansys.shell7gen
                    delete(shell7_name);  % shell7 not requested, delete temp file
                end
                cd(old_dir);
            catch ME
                cd(old_dir);
                rethrow(ME);
            end
        end
    end
end

function cb_gen_plot3d(cbo,~)
    app = guidata(cbo);  % load application data
    blade = NuMAD_appdata('get','blade');
    
    if isempty(app.settings.job_name)
        warndlg('Please save the model before generating output.','Operation Not Permitted');
        return;
    end
    
    [~,name,~] = fileparts(app.settings.job_name);
    filename = fullfile(app.settings.job_path,[name '.p3d']);
    %jcb: should we instead ask user for filename?
    
    % check that all newspanloc values are within the blade span
    spanloc = [app.station.LocationZ];
    
    values = app.plot3d.newspanloc;
    outofrange = (values<spanloc(1)) | (values>spanloc(end));
    for k=1:numel(values)
        if ~isempty(find(values(k)==spanloc,1))
            outofrange(k) = 1;
        end
    end
    if any(outofrange)
        message = ['Cannot interpolate at the following spans: ', sprintf('%g ',values(outofrange))];
        warndlg(message,'Plot3D output warning');
    end
    values(outofrange) = [];
    app.plot3d.newspanloc = sort(values);
    
    write_plot3d(app,blade,filename);
end

function cb_runbpe(cbo,~)
    app = guidata(cbo);  % load application data
    data.station = app.station;
    data.BladeRotation = app.BladeRotation;
    
    previous_dir = pwd;
    try
        cd(app.settings.job_path);
        NuMAD2BPE2FASTBlade(cbo,data,app.settings);
        
        %     delete FASTBlade.dat
        delete utlsuite*
        delete BPE_SectionData.mat
        delete bmodes.out
        delete shell7bpe.src
        %     delete displacement.txt
        delete output.txt
        disp('BPE-related files have been deleted.')
        cd(previous_dir);
    catch ME
        cd(previous_dir);
        rethrow(ME);
    end
end

function cb_runprecomp(cbo,~)
    app = guidata(cbo);  % load application data
    app.blade = NuMAD_appdata('get','blade'); %brr
    data.station = app.station;
    data.shearweb = app.shearweb;
    data.active.list = app.complist;
    data.BladeRotation = app.BladeRotation;
    data.matdb = app.matdb;
    data.blade = app.blade;
    if app.batchrun
        data.batchmodelist = app.batchmodelist;
    end
    
    for i=1:numel(data.station)
        if ~any(data.station(i).dp==0)
            errordlg('Each station must contain a material D.P. at the leading edge in order for NuMAD2PreComp to work.','Error: NuMAD2PreComp');
            return;
        end
    end
    
    previous_dir = pwd;
    try
        cd(app.settings.job_path);
        bmodesFrequencies = NuMAD2PreComp2FASTBlade(cbo,data,app.settings,app.batchrun); %#ok<NASGU>
        if app.batchrun
            save('bmodesFrequencies','bmodesFrequencies');
        end
        cd(previous_dir);
    catch ME
        cd(previous_dir);
        rethrow(ME);
    end
    
end

function cb_xls2nmd(cbo,~)
    app = guidata(cbo);  % load application data
    
    if isempty(app.settings.job_path)
        if ispc
            starting_dir = getenv('USERPROFILE');
        elseif isunix
            starting_dir = getenv('HOME');
        else
            starting_dir = '';
        end
    else
        starting_dir = app.settings.job_path;
    end

    [fn,path2xls] = uigetfile( ...
    {'*.xls;*.xlsx', 'Excel Files (*.xls,*.xlsx)';...
     'NuMAD.xls;NuMAD.xlsx', 'NuMAD.xls(x)';...
     '*.*', 'All files (*.*)'},...
     'XLS2NMD: Select NuMAD.xls file',starting_dir);

    if isequal(fn,0)
        % the user canceled file selection
        return
    end
    
    xlsname = fullfile(path2xls,fn);

    previous_dir = pwd;
    try
        cd(path2xls);
        output_file = xls2nmd(xlsname);
        cd(previous_dir);
        helpdlg('xls2nmd completed without error. Generated model will now open.','XLS-to-NuMAD');
    catch ME
        cd(previous_dir);
        warndlg('xls2nmd encountered an error','XLS-to-NuMAD');
        rethrow(ME);
    end
    % open the model
    openmodel(path2xls,output_file,cbo);
end

function cb_checkbladedata(cbo,~)
    app = guidata(cbo);  % load application data
    
    TotalStations = numel(app.station);
    TotalShearwebs = numel(app.shearweb);
    SkinAreas(1:TotalStations-1) = struct('startIB',[],'endIB',[],'startOB',[],'endOB',[],'Material',{''});
    matlist = {};
    for kStation = 1:numel(SkinAreas)
        stationIB = app.station(kStation);
        for kdp = 1:numel(stationIB.dp)-1
            switch stationIB.dptype{kdp}
                case {'single','double'}
                    % single and double are equivalent on the area inboard edge
                    % start and end are current and next DP
                    SkinAreas(kStation).startIB(end+1) = kdp;
                    SkinAreas(kStation).endIB(end+1)   = kdp+1;
                    [SkinAreas(kStation).Material{end+1} matlist{end+1}] = deal(stationIB.sm{kdp});
                case {'flare','hourglass'}
                    % flare and hourglass are equivalent on the area inboard edge
                    % start and end of first area is current DP
                    % start and end of next area is current and next DP
                    SkinAreas(kStation).startIB(end+(1:2)) = [kdp kdp];
                    SkinAreas(kStation).endIB(end+(1:2))   = [kdp kdp+1];
                    [SkinAreas(kStation).Material{end+1} matlist{end+1}] = deal(stationIB.dpmaterial{kdp});  %jcb:  append (end+1) dp material
                    [SkinAreas(kStation).Material{end+1} matlist{end+1}] = deal(stationIB.sm{kdp});  %jcb:  append (end+1) segment material
            end
        end
        app.matlist = unique(matlist);
        
        stationOB = app.station(kStation+1);
        for kdp = 1:numel(stationOB.dp)-1
            switch stationOB.dptype{kdp}
                case {'single','flare'}
                    % single and flare are equivalent on the area outboard edge
                    % start and end are current and next DP
                    SkinAreas(kStation).startOB(end+1) = kdp;
                    SkinAreas(kStation).endOB(end+1)   = kdp+1;
                case {'double','hourglass'}
                    % double and hourglass are equivalent on the area outboard edge
                    % start and end of first area is current DP
                    % start and end of next area is current and next DP
                    SkinAreas(kStation).startOB(end+(1:2)) = [kdp kdp];
                    SkinAreas(kStation).endOB(end+(1:2))   = [kdp kdp+1];
            end
        end
    end

    app.SkinAreas = SkinAreas;
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
    
    app = draw_skinareas(app);
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end
    
    check_results = '';
    for kStation = 1:numel(SkinAreas)
        nSegmentsIB = numel(SkinAreas(kStation).startIB);
        nSegmentsOB = numel(SkinAreas(kStation).startOB);
        if nSegmentsIB ~= nSegmentsOB
            check_results = sprintf('%sDP mismatch: Station %d = %d segments,',check_results,kStation,nSegmentsIB);
            check_results = sprintf('%s Station %d = %d segments\n',check_results,kStation+1,nSegmentsOB);
        end
        if any(strcmp('**UNSPECIFIED**',SkinAreas(kStation).Material))
            check_results = sprintf('%sUnspecified material(s): Station %d.\n',check_results,kStation);
        end
        if isequal(app.station(kStation).TEtype,'flat') ...
                && (app.station(kStation).dp(2)~=-1 || app.station(kStation).dp(end-1)~=1)
            check_results = sprintf('%sStation %d: a flatback airfoil must have DPs at both corners of the flat.\n',check_results,kStation);  
        end
    end
    for ksw = 1:TotalShearwebs
        if strcmp('**UNSPECIFIED**',app.shearweb(ksw).Material)
            check_results = sprintf('%sUnspecified material: Station %d shear web.\n',check_results,app.shearweb(ksw).BeginStation);
        end
    end
    if isempty(check_results)
        if ~app.modifysmd && ~app.batchrun
            helpdlg('Model data verified.','Model Check');
        end
        app.checkpassed = true;
    else
        if ~app.batchrun
            warndlg(check_results,'Model Check');
        else
            fprintf(check_results);
        end
        app.checkpassed = false;
    end
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end 
end

function cb_checkmatdata(cbo,~)
    app = guidata(cbo);  % load application data
    
    %Determine which composite materials are used in the model
    if isempty(app.shearweb)
        compsInModel = {};
    else
        compsInModel = {app.shearweb.Material}';
    end
    for k = 1:numel(app.SkinAreas)
        compsInModel = [compsInModel; transpose(app.SkinAreas(k).Material)];  
    end
    compsInModel = unique(compsInModel);
    
    matlist = {app.matdb.name}';
%     mattype = {app.matdb.type}';
%     isotropic =  strcmp('isotropic',mattype);
%     orthotropic = strcmp('orthotropic',mattype);
%     composite = strcmp('composite',mattype);
    
    % Determine which isotropic and orthotropic materials are used in the model
    isoorthoInModel = {};
    for kcomp = 1:numel(compsInModel)
        n = strcmp(compsInModel(kcomp),matlist);
        if ~any(n)
            errordlg(sprintf('Material "%s" not found in database.',compsInModel{kcomp}),'Error');
            error('Material "%s" not found in database.',compsInModel{kcomp});
        end
        mat = app.matdb(n);
        layerNames = {mat.layer.layerName}';
        isoorthoInModel = unique([isoorthoInModel; layerNames]);
    end
    
    fcfields = {'xten','xcmp','yten','ycmp','zten','zcmp','xy','yz','xz',...
        'xycp','yzcp','xzcp','xzit','xzic','yzit','yzic','g1g2',...
        'etal','etat','alp0'};
    check_results = {};
    missingFCProps = zeros(numel(isoorthoInModel),length(app.ansys.FailureCriteria));
    fcselections = cell2mat(app.ansys.FailureCriteria(:,2));
    for kmat = 1:numel(isoorthoInModel)
        fcmiss = '';
        n = strcmp(isoorthoInModel{kmat},matlist);
        mat = app.matdb(n);
        for kfc = 1:numel(fcfields)
            fcprop = fcfields{kfc};
            fcvalues{kfc} = mat.(fcprop);
        end
        fcempty = cellfun('isempty',fcvalues);
        if any(fcempty(1:9))  % Max S
            missingFCProps(kmat,1:2) = 1; 
            if any(fcselections(1:2))
                fcmiss = [fcmiss 'ESMax, '];
            end
        end  
        if any(fcempty(1:12))  % Tsai-Wu
            missingFCProps(kmat,3:4) = 1;
            if any(fcselections(3:4))
                fcmiss = [fcmiss 'Tsai-Wu, '];
            end
        end
        if any(fcempty([1:4,6:7]))  % Hashin
            missingFCProps(kmat,5:6) = 1; 
            if any(fcselections(5:6))
                fcmiss = [fcmiss 'Hashin, '];
            end
        end
        if any(fcempty([1:4,7,13:16]))  % Puck
            missingFCProps(kmat,7:8) = 1;
            if any(fcselections(7:8))
                fcmiss = [fcmiss 'Puck, '];
            end
        end
        if any(fcempty([1:4,7,17:20]))  % LaRc
            missingFCProps(kmat,9:12) = 1; 
            if any(fcselections(9:12))
                fcmiss = [fcmiss 'LaRc, '];
            end
        end
        if any(fcempty(1:20))  % User
            missingFCProps(kmat,13:21) = 1; 
            if any(fcselections(13:21))
                fcmiss = [fcmiss 'User, '];
            end
        end
        
        if ~isempty(fcmiss)
            fcmiss(end-1:end) = [];  %remove final comma
            msg = sprintf('Material "%s" lacks properties needed for failure criteria %s.',mat.name,fcmiss);
            %check_results = sprintf('%s%s\n',check_results,msg);
            check_results{end+1} = msg;
        end

    end
    
    app.missingFCProps = missingFCProps;
    
    if isempty(check_results)
        
    else
        h = figure('Name','Material Check',...
            'NumberTitle','off',...
            'NextPlot','new',...
            'Menubar','none',...
            'Resize','off',...
            'Position',[10000 10000 600 200]);
        uicontrol(h,...
            'Style','listbox',...
            'Max',2,...  % setting Max=2 Min=0 allows multiply lines of text
            'Min',0,...
            'HorizontalAlignment','left',...
            'String',check_results,...
            'BackgroundColor','white',...
            'Position',[10 10 580 180]);
        movegui(h,'center');
    end
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app',app);
    end 
end

function cb_help(~,~,option)
    switch option
        case 'userguide'
            try
                open('NuMAD_UserGuide_RandA.pdf')
            catch %#ok<CTCH>
                fprintf('\nCould not open user guide.\n');
                warndlg('Could not open user guide.','Help');
            end
        case 'userguideweb'
            stat = web('http://energy.sandia.gov/?tag=numad','-browser');
            switch stat
                case 1
                    fprintf('\nBrowser could not be found.\n')
                    warndlg('Browser could not be found.','Help');
                case 2
                    fprintf('\nBrowser found, but could not be launched.\n')
                    warndlg('Browser found, but could not be launched.','Help');
            end
        case 'aboutnumad'
            h = figure('Name','About NuMAD',...
                'BackingStore','off',...
                'DockControls','off',...
                'HandleVisibility','callback',...
                'Integerhandle','off',...
                'Menubar','none',...
                'NumberTitle','off',...
                'Resize','off',...
                'Visible','on',...
                'NextPlot','new',...
                'ButtonDownFcn','close(gcbf)');
            about_numad = {'';
                '**********************************************';
                '*       Part of the SNL NuMAD Toolbox        *';
                '* Developed by Sandia National Laboratories, *'
                '*           Wind Energy Technologies         *';
                '*                                            *';
                '* See license.txt for disclaimer information *';
                '**********************************************'};
            g=uicontrol(h,'Style','text',...
                'max',100,...
                'min',1,...
                'string',about_numad,...
                'Units','normalized',...
                'Position',[0.05 0.05 0.9 0.9],...
                'FontName','Courier New',...
                'HorizontalAlignment','center',...
                'ButtonDownFcn','close(gcbf)');

    end
end

function pointer = getPointerShapes()
% manual process:
%   #load 8-bit rgb image#
%   m=mean(img,3);  % average of RGB channels
%   p=nan(16); p(m==0)=1; p(m==255)=2;  % black=1, white=2, otherwise=NaN

% mouse pointer shape for zoom
pointer.zoom = [NaN,NaN,NaN,NaN,1,1,1,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,1,1,2,2,2,2,1,1,NaN,NaN,NaN,NaN,NaN,NaN;NaN,1,2,2,2,NaN,NaN,2,2,2,1,NaN,NaN,NaN,NaN,NaN;NaN,1,2,NaN,NaN,NaN,NaN,NaN,NaN,2,1,NaN,NaN,NaN,NaN,NaN;1,2,2,NaN,NaN,NaN,NaN,NaN,NaN,2,2,1,NaN,NaN,NaN,NaN;1,2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,2,1,NaN,NaN,NaN,NaN;1,2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,2,1,NaN,NaN,NaN,NaN;1,2,2,NaN,NaN,NaN,NaN,NaN,NaN,2,2,1,NaN,NaN,NaN,NaN;NaN,1,2,NaN,NaN,NaN,NaN,NaN,NaN,2,1,NaN,NaN,NaN,NaN,NaN;NaN,1,2,2,2,NaN,NaN,2,2,2,1,NaN,NaN,NaN,NaN,NaN;NaN,NaN,1,1,2,2,2,2,1,1,2,1,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,1,1,1,1,NaN,NaN,1,1,1,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,1,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,1,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,1;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1;];
% mouse pointer shape for rotate 3d (orbit)
pointer.rotate = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,2,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,2,1,2,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,2,1,2,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,2,1,2,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,2,2,2,2,2,2,2,NaN,NaN,NaN,NaN;NaN,NaN,NaN,2,2,1,1,1,1,1,1,1,2,2,NaN,NaN;NaN,NaN,2,1,1,1,2,2,2,2,2,1,1,1,2,NaN;NaN,2,1,1,2,2,NaN,2,1,2,NaN,2,2,1,1,2;NaN,2,1,2,NaN,NaN,NaN,2,1,2,NaN,NaN,NaN,2,1,2;NaN,2,1,2,NaN,NaN,NaN,2,1,2,NaN,NaN,NaN,2,1,2;NaN,2,1,1,2,2,NaN,NaN,2,NaN,NaN,2,2,1,1,2;NaN,NaN,2,1,2,1,2,NaN,2,NaN,2,1,2,1,2,NaN;NaN,NaN,NaN,2,1,1,2,2,1,2,2,1,1,2,NaN,NaN;NaN,NaN,2,1,1,1,2,2,1,2,2,1,1,1,2,NaN;NaN,NaN,NaN,2,2,2,NaN,2,1,2,NaN,2,2,2,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,2,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
end