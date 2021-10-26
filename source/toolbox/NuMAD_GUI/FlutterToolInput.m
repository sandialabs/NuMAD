function FlutterToolInput(varargin)
%FlutterToolInput  User interface to flutter tool
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   FlutterToolInput(cbo), in standard usage, is used as a menu callback 
%   within NuMAD_main which launches the flutter tool user interface. 
%   Usage: uimenu('label','Flutter Tool','callback',@FlutterToolInput)
%
%   FlutterToolInput(FILENAME) opens flutter tool interface outside NuMAD
%      FILENAME  = filename of a NuMAD project file

%===== INITIALIZATION =====================================================

if ischar(varargin{1})
    % assume filename argument
    nmdfn = varargin{1};
    [app.station, shearweb, active, ansys, BladeRotation, blade, app.flutterinput] = readNuMADinput(nmdfn);
    mp = get(0,'MonitorPositions');
    xpos = (mp(1,3)-800)/2;
    ypos = (mp(1,4)-300)/2;
    app.xy = [xpos ypos];
    app.job_path = fileparts(nmdfn);
else
    cbo = varargin{1};
    callapp = guidata(cbo);
    if isempty(callapp.settings.job_name)
        warndlg('Please save the model before running this tool.','Operation Not Permitted');
        return;
    end
    app.caller = callapp.fh;  % store figure handle provided by calling script
    app.numadpath = callapp.numadpath;
    app.userpath = callapp.userpath;
    app.job_name = callapp.settings.job_name;
    app.job_path = callapp.settings.job_path;
    app.station = callapp.station;
    if isfield(callapp,'flutterInput')
        app.flutterInput = callapp.flutterInput;
    else
        app.flutterInput.fstFile = '';
        app.flutterInput.bladeFile = '';
        app.flutterInput.aeroFile = '';
        app.flutterInput.outFile = 'name.out';
        app.flutterInput.LCS = [];
        app.flutterInput.OmegaStart = [];
        app.flutterInput.OmegaIncrement = [];
        app.flutterInput.OmegaEnd = [];
    end
    app.settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
    app.xy = app.settings.xy_flutterinput;
end

app.debugging = 0;  % with debugging==1, the data structure is sent to the
% workspace after the gui is constructed
%===== GUI CONSTRUCTION ===================================================
c=cgui([],'init');  c.fig.xy = app.xy;
c=cgui(c,'figure','Flutter Tool',[630 470],'guiFlutterTool');
c=cgui(c,'panel','',[10 10 610 450]);
paneltop = c.ctrl.vpos;
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 80 500]);
c=cgui(c,'input','FAST file:',app.flutterInput.fstFile,'fstFile',@cb_filenameinput);
c=cgui(c,'input','FAST blade:',app.flutterInput.bladeFile,'bladeFile',@cb_filenameinput);
c=cgui(c,'input','AeroDyn file:',app.flutterInput.aeroFile,'aeroFile',@cb_filenameinput);
c=cgui(c,'input','Results file:',app.flutterInput.outFile,'outFile',@cb_filenameinput);
c=cgui(c,'colEdges','',[10 120]);
c=cgui(c,'label','Lift curve slope:','','','');
c=cgui(c,'vpos','value',paneltop);
c=cgui(c,'colEdges','',[510 600]);
c=cgui(c,'push','Browse','','',@cb_browsefastfile);
c=cgui(c,'push','Browse','','',@cb_browsebladefile);
c=cgui(c,'push','Browse','','',@cb_browseaerofile);
c=cgui(c,'push','Browse','','',@cb_browseoutfile);

tableColumns = {'Airfoil files', 'LCS'};
tableData = {};
[pn,fn,ext] = fileparts(app.flutterInput.aeroFile);
if isempty(pn)
    aeroFile = fullfile(app.job_path,app.flutterInput.aeroFile);
else
    aeroFile = app.flutterInput.aeroFile;
end
if isequal(2,exist(aeroFile,'file'))
    try
        app.ad = readFastAD(aeroFile);
        tableData = app.ad.FoilNm;
    catch
        warndlg('Could not read AeroDyn file.','Flutter Tool');
    end
end
LCSn = numel(app.flutterInput.LCS);
if LCSn>0
    if isequal(LCSn,size(tableData,1))
        tableData(:,2) = num2cell(app.flutterInput.LCS);
    else
        % size of LCS array does not match
        warndlg('Size of stored LCS array does not match number of airfoils in AeroDyn file. Stored LCS values discarded.','Flutter Tool');
        if ~isempty(tableData)
            tableData(:,2) = {6.28};  % 2*pi = 6.28
        end
    end
end
app.th=uitable('Parent',c.panel.h,'Position',[80 (paneltop-324) 410 200],...
    'ColumnWidth',{255 55},'ColumnEditable',[false true],...
    'Data',tableData,'ColumnName',tableColumns); 

c=cgui(c,'vpos','value',paneltop-150);
c=cgui(c,'rowHeight','',60);
c=cgui(c,'colEdges','',[510 600]);
c=cgui(c,'push','Auto Populate','','',@cb_autolcs);

c=cgui(c,'vpos','value',paneltop-340);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 100 200]);
c=cgui(c,'input','RPM Start',app.flutterInput.OmegaStart,'OmegaStart',@cb_num);
c=cgui(c,'input','RPM End',app.flutterInput.OmegaEnd,'OmegaEnd',@cb_num);
c=cgui(c,'input','RPM Increment',app.flutterInput.OmegaIncrement,'OmegaIncrement',@cb_num);

c=cgui(c,'vpos','value',paneltop-360);
c=cgui(c,'rowHeight','',40);
c=cgui(c,'colEdges','',[300 490]);
c=cgui(c,'push','GO','','',@cb_goflutter);

app.fh = c.fig.h;
set(app.fh,'CloseRequestFcn',@cb_save);
guidata(app.fh,app);  % store the application data as guidata

if app.debugging
    assignin('base','app_fluttertool',app); 
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

function cb_filenameinput(cbo,~)
    app = guidata(cbo);  % load application data
    tag = get(cbo,'Tag');
    filenameinput = get(cbo,'String');
    if isempty(filenameinput)
        filenameinput = '';
    end
    [pn,fn,ext] = fileparts(filenameinput);
    if isequal('',pn) && ~isequal('',fn)
        filenametest = fullfile(app.job_path,[fn,ext]);
    else
        filenametest = filenameinput;
    end
    
    if isequal(tag,'outFile')
        app.flutterInput.(tag) = filenameinput;
    else
        if isempty(filenameinput) || isequal(2,exist(filenametest,'file'))  % 2 indicates name is a file
            app.flutterInput.(tag) = filenameinput;
            set(cbo,'BackgroundColor','white');
        else
            set(cbo,'BackgroundColor',[1 1 0.8]);
        end
    end
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_fluttertool',app);
    end
end

function cb_browsefastfile(cbo,~)
    app = guidata(cbo);  % load application data
    [fn pn] = uigetfile( ...
        {'*.*', 'All Files (*.*)';'*.fst','FAST file (*.fst)'}, ...
        'Select FAST input file',app.job_path);
    if isequal(fn,0) || isequal(pn,0)
        return
    else
        if isequal(pn,app.job_path)
            app.flutterInput.fstFile = fn;
        else
            app.flutterInput.fstFile = fullfile(pn,fn);
        end
        guidata(app.fh,app);  % store the updated guidata
        if app.debugging
            assignin('base','app_fluttertool',app);
        end
        uictrl = get_uictrl(gcbf,'fstFile');
        set(uictrl.h,'String',app.flutterInput.fstFile);
    end
end

function cb_browsebladefile(cbo,~)
    app = guidata(cbo);  % load application data
    [fn pn] = uigetfile( ...
        {'*.*', 'All Files (*.*)';'*.dat','FAST blade file (*.dat)'}, ...
        'Select FAST blade file',app.job_path);
    if isequal(fn,0) || isequal(pn,0)
        return
    else
        if isequal(pn,app.job_path)
            app.flutterInput.bladeFile = fn;
        else
            app.flutterInput.bladeFile = fullfile(pn,fn);
        end
        guidata(app.fh,app);  % store the updated guidata
        if app.debugging
            assignin('base','app_fluttertool',app);
        end
        uictrl = get_uictrl(gcbf,'bladeFile');
        set(uictrl.h,'String',app.flutterInput.bladeFile);
    end
end

function cb_browseaerofile(cbo,~)
    app = guidata(cbo);  % load application data
    [fn pn] = uigetfile( ...
        {'*.*', 'All Files (*.*)';'*.ipt','AeroDyn file (*.ipt)'}, ...
        'Select AeroDyn file',app.job_path);
    if isequal(fn,0) || isequal(pn,0)
        return
    else
        if isequal(pn,app.job_path)
            app.flutterInput.aeroFile = fn;
        else
            app.flutterInput.aeroFile = fullfile(pn,fn);
        end
        try
            app.ad = readFastAD(fullfile(pn,fn));
            tableData = app.ad.FoilNm;
        catch
            warndlg('Could not read AeroDyn file.','Flutter Tool');
            tableData = {};
        end
        guidata(app.fh,app);  % store the updated guidata
        if app.debugging
            assignin('base','app_fluttertool',app);
        end
        % update textbox
        uictrl = get_uictrl(gcbf,'aeroFile');
        set(uictrl.h,'String',app.flutterInput.aeroFile);
        % update table
        originalTableData = get(app.th,'Data');
        LCSn = size(originalTableData,1);
        if isequal(LCSn,size(tableData,1))
            tableData(:,2) = originalTableData(:,2);
        else
            % size of LCS array does not match
            warndlg('Size of stored LCS array does not match number of airfoils in AeroDyn file. Stored LCS values discarded.','Flutter Tool');
            if ~isempty(tableData)
                tableData(:,2) = {6.28};  % 2*pi = 6.28
            end
        end
        set(app.th,'Data',tableData);
    end
end

function cb_browseoutfile(cbo,~)
    app = guidata(cbo);  % load application data
    [fn pn] = uiputfile( ...
        {'*.*', 'All Files (*.*)';'*.out','Results file (*.out)'}, ...
        'Select output file for results',fullfile(app.job_path,'name.out'));
    if isequal(fn,0) || isequal(pn,0)
        return
    else
        if isequal(pn,app.job_path)
            app.flutterInput.outFile = fn;
        else
            app.flutterInput.outFile = fullfile(pn,fn);
        end
        guidata(app.fh,app);  % store the updated guidata
        if app.debugging
            assignin('base','app_fluttertool',app);
        end
        uictrl = get_uictrl(gcbf,'outFile');
        set(uictrl.h,'String',app.flutterInput.outFile);
    end
end

function cb_autolcs(cbo,~)
    app = guidata(cbo);     % retrieve application data
    
    [pn,fn,ext] = fileparts(app.flutterInput.aeroFile);
    if isempty(pn);
        pn = app.job_path;
    end
    aeroFile = fullfile(pn,[fn,ext]);
    try
        app.ad = readFastAD(aeroFile);
    catch
        warndlg('Could not read AeroDyn file.','Flutter Tool');
        return
    end
    
    LCS = zeros(app.ad.NumFoil,1);
    for kf = 1:length(LCS)
        foilname = app.ad.FoilNm{kf};
        foilname = strrep(foilname,'"','');
        try
            foildata = readAirfoilData(fullfile(pn,foilname));
        catch
            fprintf('Could not read %s.\n',fullfile(pn,foilname));
            continue
        end
        AoA = foildata.AoA;
        CL = foildata.CL(1,:); % hard-coded to first table in multi-table files
        
        LCS(kf) = findLCS(AoA,CL);
    end
    
    app.flutterInput.LCS = LCS;
    tableData = app.ad.FoilNm;
    tableData(:,2) = num2cell(LCS);
    set(app.th,'Data',tableData);
end

function LCS = findLCS(AoA,CL)
    k0 = find(AoA>=0,1);
    kf = find(AoA>=12,1);
    threshold = 0.02;
    LCS = nan;
    for k=(k0+2):kf
        [p,s] = polyfit(AoA(k0:k),CL(k0:k),1);
        if s.normr < threshold
            LCS = p(1)*180/pi;
        end
        if (0)
            CLfit = polyval(p,AoA(k0:k));
            plot(AoA(k0:k),CL(k0:k),AoA(k0:k),CLfit);
            title(sprintf('slope=%g, normr=%g',p(1)*180/pi,s.normr));
            pause(0.5);
        end
    end
end

function cb_num(cbo,~)
% this callback checks numeric input and updates app data
    app = guidata(cbo);     % retrieve application data
    tag = get(cbo,'tag');  % get the object tag (data struct fieldname)
    % get the stored data associated with the fieldname
    val = app.flutterInput.(tag);  % use dynamic fieldname to get value
    str = get(cbo,'string');     % get the new input
    num = str2double(str);  % try to convert to number
    if isnan(num)
        % input was not valid => revert to stored value
        str = sprintf('%g',val);
        set(cbo,'string',str);
    else
        % input valid => save new value
        app.flutterInput.(tag) = num;
        %blade.(tag)(k) = app.station(k).(tag);
        guidata(app.fh,app);  % store the updated guidata
        if app.debugging
            assignin('base','app_fluttertool',app);
        end
    end
end

function cb_goflutter(cbo,~)
    app = guidata(cbo);  % load application data
    
    flutterInput.fstFile = app.flutterInput.fstFile;
    flutterInput.bladeFile = app.flutterInput.bladeFile;
    flutterInput.aeroFile = app.flutterInput.aeroFile;
    flutterInput.outFile = app.flutterInput.outFile;
    
    tableData = get(app.th,'Data');
    app.flutterInput.LCS = cell2mat(tableData(:,2));
    flutterInput.LCS = app.flutterInput.LCS;
    
    callapp = guidata(app.caller);
    flutterInput.pitchAxisDomain = [callapp.station.LocationZ];
    flutterInput.pitchAxisVal = [callapp.station.Xoffset];
    len = length(flutterInput.pitchAxisVal);
    flutterInput.numadBladeLen = (callapp.station(len).LocationZ-callapp.station(1).LocationZ);
    
    OmegaStart = app.flutterInput.OmegaStart;
    OmegaEnd = app.flutterInput.OmegaEnd;
    OmegaIncrement = app.flutterInput.OmegaIncrement;
    flutterInput.OmegaArray = OmegaStart:OmegaIncrement:OmegaEnd; 
    
    savefile = fullfile(app.job_path,'flutterInput.mat');
    save(savefile,'flutterInput');
    
    tempPath = pwd;
    cd(app.job_path);
    try
        [convergedFreq,convergedDamp]=feaAutoAllModes(flutterInput,1.0e-3,'F');
    catch ME
        cd(tempPath);
        rethrow(ME);
    end
    cd(tempPath);
end

function cb_save(cbo,~)
    app = guidata(cbo);  % load application data
    
    tableData = get(app.th,'Data');
    if size(tableData,1)>0
        app.flutterInput.LCS = cell2mat(tableData(:,2));
    else
        app.flutterInput.LCS = [];
    end
    
    if isobject(app.caller) % verify main NuMAD window is still open
        callapp = guidata(app.caller);
        callapp.flutterInput = app.flutterInput;
        guidata(callapp.fh,callapp);  % store the updated application data
        if callapp.debugging
            assignin('base','app',callapp);
        end
    end
    
    cb_close(cbo,[]);
end

function cb_close(cbo,~)
    try
        app = guidata(cbo);  % load application data
        position = get(app.fh,'Position');
        settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
        if ~isequal(settings.xy_flutterinput,position(1:2))
            settings.xy_flutterinput = position(1:2);
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
