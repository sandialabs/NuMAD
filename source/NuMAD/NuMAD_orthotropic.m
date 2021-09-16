function NuMAD_orthotropic(mat,caller,cb_matsave)
%NUMAD_ORTHOTROPIC  User interface for orthotropic material definition 
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_orthotropic launches the interface which allows the user to 
%   define the properties of orthotropic materials.
%

%===== INITIALIZATION =====================================================
if nargin==3
    if ~isstruct(mat) || ~isfield(mat,'type') || ~strcmp(mat.type,'orthotropic')
        errordlg('Argument of this function must be an orthotropic material data structure.','Programming Error');
        return;
    end
else
    errordlg('Script must be called from main materials window.','Programming Error');
    return;
end

callapp = guidata(caller);
app.userpath = callapp.userpath;

if strcmp(mat.name,'__NewMatFlag')
    orgname = mat.name;  % remember the original material name in case the user changes it
    mat.name = 'new orthotropic';
    mat.dens = 1000;
    mat.ex = 1e9;
    mat.ey = 1e9;
    mat.ez = 1e9;
    mat.gxy = 1e9;
    mat.gyz = 1e9;
    mat.gxz = 1e9;
    mat.prxy = 0.3;
    mat.pryz = 0.3;
    mat.prxz = 0.3;
else
    orgname = mat.name;  % remember the original material name in case the user changes it
end

app.debugging = 0;  % with debugging==1, the data structure is sent to the 
                    % workspace after the gui is constructed
%===== GUI CONSTRUCTION ===================================================
units = unitconst();

app.settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
c=cgui([],'init');  c.fig.xy = app.settings.xy_orthotropic;
c=cgui(c,'figure','NuMAD - modify orthotropic',[500 700],'guiNuMAD-orthotropic');
set(c.fig.h,'CloseRequestFcn',@cb_close);

uicontrol(c.fig.h,...
            'Style','popupmenu',...
            'String',{'Elastic properties','Failure properties'},...
            'Value',1,...
            'Tag','PropertySet',...
            'BackgroundColor','white',...
            'Position',[10 660 120 22],...
            'Callback',@cb_propertyset);
app.crit(1) = uicontrol(c.fig.h,...
            'Style','text',...
            'String','ESMax',...
            'Tag','esmax',...
            'Position',[150 660 50 20]);
app.crit(2) = uicontrol(c.fig.h,...
            'Style','text',...
            'String','Tsai-Wu',...
            'Tag','tsaiwu',...
            'Position',[200 660 50 20]);
app.crit(3) = uicontrol(c.fig.h,...
            'Style','text',...
            'String','Puck',...
            'Tag','puck',...
            'Position',[250 660 50 20]);
app.crit(4) = uicontrol(c.fig.h,...
            'Style','text',...
            'String','Hashin',...
            'Tag','hashin',...
            'Position',[300 660 50 20]);
app.crit(5) = uicontrol(c.fig.h,...
            'Style','text',...
            'String','LaRc',...
            'Tag','larc',...
            'Position',[350 660 50 20]);
app.crit(6) = uicontrol(c.fig.h,...
            'Style','text',...
            'String','User',...
            'Tag','user',...
            'Position',[400 660 50 20]);
set(app.crit,'Visible','off');
uicontrol(c.fig.h,...
            'Style','pushbutton',...
            'String','Save Changes',...
            'Position',[10 10 320 22],...
            'Callback',{@cb_save,cb_matsave});
uicontrol(c.fig.h,...
            'Style','pushbutton',...
            'String','Discard Changes',...
            'Position',[340 10 150 22],...
            'Callback',@cb_close);   
% jcb: rather that just using 'close', we should maybe check if anything 
% has changed and confirm with user that they wish to abandon changes

c=cgui(c,'panel','Elastic properties',[10 40 480 600]);
app.panel(1) = c.panel.h;
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 160     360]);
c=cgui(c,'input','Material Name:',mat.name,'name',@cb_text);
c=cgui(c,'colEdges','',[10 160 270 360]);
c=cgui(c,'dnum','Young''s Modulus (Ex):',mat.ex,'ex',@cb_dnum,units.stress);
c=cgui(c,'dnum','Young''s Modulus (Ey):',mat.ey,'ey',@cb_dnum,units.stress);
c=cgui(c,'dnum','Young''s Modulus (Ez):',mat.ez,'ez',@cb_dnum,units.stress);
c=cgui(c,'dnum','Shear Modulus (Gxy):',mat.gxy,'gxy',@cb_dnum,units.stress);
c=cgui(c,'dnum','Shear Modulus (Gyz):',mat.gyz,'gyz',@cb_dnum,units.stress);
c=cgui(c,'dnum','Shear Modulus (Gxz):',mat.gxz,'gxz',@cb_dnum,units.stress);
c=cgui(c,'num' ,'Major Poisson''s Ratio (XY):',mat.prxy,'prxy',@cb_num);
c=cgui(c,'num' ,'Major Poisson''s Ratio (YZ):',mat.pryz,'pryz',@cb_num);
c=cgui(c,'num' ,'Major Poisson''s Ratio (XZ):',mat.prxz,'prxz',@cb_num);
c=cgui(c,'dnum','Mass Density:',mat.dens,'dens',@cb_dnum,units.density);
c=cgui(c,'rowHeight','',100);
c=cgui(c,'colEdges','',[10 80 470]);
c=cgui(c,'textbox','Reference:',mat.reference,'reference',@cb_multilinetext);

c=cgui(c,'panel','Failure properties',[10 40 480 600]);
app.panel(2) = c.panel.h;
set(c.panel.h,'Visible','off');
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 230 340 390]);
c=cgui(c,'dnum','Allowable X tensile stress (XTEN)',mat.xten,'xten',@cb_dnum2,units.stress);
c=cgui(c,'dnum','Allowable X compressive stress (XCMP)',mat.xcmp,'xcmp',@cb_dnum2,units.stress);
c=cgui(c,'dnum','Allowable Y tensile stress (YTEN)',mat.yten,'yten',@cb_dnum2,units.stress);
c=cgui(c,'dnum','Allowable Y compressive stress (YCMP)',mat.ycmp,'ycmp',@cb_dnum2,units.stress);
c=cgui(c,'dnum','Allowable Z tensile stress (ZTEN)',mat.zten,'zten',@cb_dnum2,units.stress);
c=cgui(c,'dnum','Allowable Z compressive stress (ZCMP)',mat.zcmp,'zcmp',@cb_dnum2,units.stress);
c=cgui(c,'dnum','Allowable XY shear stress (XY)',mat.xy,'xy',@cb_dnum2,units.stress);
c=cgui(c,'dnum','Allowable YZ shear stress (YZ)',mat.yz,'yz',@cb_dnum2,units.stress);
c=cgui(c,'dnum','Allowable XZ shear stress (XZ)',mat.xz,'xz',@cb_dnum2,units.stress);
c=cgui(c,'num','Tsai-Wu XY coupling coeff. (XYCP)',mat.xycp,'xycp',@cb_num2);
c=cgui(c,'num','Tsai-Wu YZ coupling coeff. (YZCP)',mat.yzcp,'yzcp',@cb_num2);
c=cgui(c,'num','Tsai-Wu XZ coupling coeff. (XZCP)',mat.xzcp,'xzcp',@cb_num2);
c=cgui(c,'num','Puck XZ tensile inclination (XZIT)',mat.xzit,'xzit',@cb_num2);
c=cgui(c,'num','Puck XZ compressive inclination (XZIC)',mat.xzic,'xzic',@cb_num2);
c=cgui(c,'num','Puck YZ tensile inclination (YZIT)',mat.yzit,'yzit',@cb_num2);
c=cgui(c,'num','Puck YZ compressive inclination (YZIC)',mat.yzic,'yzic',@cb_num2);
c=cgui(c,'num','Fracture toughness ratio (G1G2)',mat.g1g2,'g1g2',@cb_num2);
c=cgui(c,'num','Longitudinal friction coefficient (ETAL)',mat.etal,'etal',@cb_num2);
c=cgui(c,'num','Transverse friction coefficient (ETAT)',mat.etat,'etat',@cb_num2);
c=cgui(c,'dnum','Fracture angle under pure transverse compression (ALP0)',mat.alp0,'alp0',@cb_dnum2,units.angle);



app.fh = c.fig.h;
app.mat = mat;
app.matname = orgname;  % remember the original material name in case the user changes it
app.caller = caller;
guidata(app.fh,app);  % store the application data as guidata

if app.debugging
    assignin('base','app_ortho',app); 
end

cb_fccheck(app.fh,[]);

end %END GUI CONSTRUCTION (main function)


%===== UTILITY & CALLBACK FUNCTIONS =======================================
function cb_text(cbo,~)
% this callback reads text input and updates app data
    app = guidata(cbo);     % retrieve guidata
    tag = get(cbo,'tag'); % get the object tag (data struct fieldname)
    str = get(cbo,'string');  % get the new input
    app.mat.(tag) = str;
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_ortho',app); 
    end
end

function cb_multilinetext(cbo,~)
% this callback reads text input and updates app data
    app = guidata(cbo);   % retrieve guidata
    tag = get(cbo,'tag'); % get the object tag (data struct fieldname)
    str = get(cbo,'string');  % get the new input
    cstr = cellstr(str);  % convert from character array to cell strings
    
    % place all the rows in a single line of text
    str = '';
    for k=1:numel(cstr)
        str = [str sprintf('%s\n',cstr{k})];
    end
    str(end) = [];  % remove the final newline
    
    app.mat.(tag) = str;
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_ortho',app); 
    end
end

function cb_num(cbo,~)
% this callback checks numeric input and updates app data
    app = guidata(cbo);     % retrieve guidata
    tag = get(cbo,'tag'); % get the object tag (data struct fieldname)
    % get the stored data associated with the fieldname
    val = app.mat.(tag);  % use dynamic fieldname to get value
    str = get(cbo,'string');  % get the new input
    num = str2double(str);  % try to convert to number
    if isnan(num)
        % input was not valid => revert to stored value
        str = sprintf('%g',val);
        set(cbo,'string',str);
    else
        % input valid => save new value
        app.mat.(tag) = num;
        guidata(app.fh,app);  % store the updated guidata
        %         str = sprintf('%g',num);
        %         set(cbo,'string',str); % diplay what was stored
        if app.debugging
            assignin('base','app_ortho',app);
        end
    end
end

function cb_num2(cbo,~)
% this callback checks numeric input and updates app data
    app = guidata(cbo);     % retrieve guidata
    tag = get(cbo,'tag'); % get the object tag (data struct fieldname)
    % get the stored data associated with the fieldname
    val = app.mat.(tag);  % use dynamic fieldname to get value
    str = get(cbo,'string');  % get the new input
    num = str2double(str);  % try to convert to number
    if isnan(num)
                % input was not valid 
                if isempty(str) % clear value if empty
                    app.mat.(tag) = [];
                    guidata(app.fh,app);  % store the updated guidata
                else % otherwise revert to stored value
                str = sprintf('%g',val);
                set(cbo,'string',str);
                end
            else
        % input valid => save new value
        app.mat.(tag) = num;
        guidata(app.fh,app);  % store the updated guidata
        %         str = sprintf('%g',num);
        %         set(cbo,'string',str); % diplay what was stored
        if app.debugging
            assignin('base','app_ortho',app);
        end
    end
    cb_fccheck(cbo,[]);
end

function cb_dnum(cbo,~)
% this callback checks dimensional numeric input and updates app data
    app = guidata(cbo);     % retrieve guidata
    style = get(cbo,'style');  % get the uicontrol style
    switch style
        case 'edit'
            % get the object tag (data struct fieldname)
            tag = get(cbo,'tag');
            % get the units
            obj = findobj(app.fh,'tag',[tag '_units']);
            units = get(obj,'UserData');
            currentunit = get(cbo,'Value');
            % get the stored data associated with the fieldname
            val = app.mat.(tag);  % use dynamic fieldname to get value
            val = val * units(currentunit);  % do unit conversion
            str = get(cbo,'string');  % get the new input
            num = str2double(str);  % try to convert to number
            if isnan(num)
                % input was not valid => revert to stored value
                str = sprintf('%g',val);
                set(cbo,'string',str);
            else
                % input valid => save new value
                app.mat.(tag) = num / units(currentunit);
                guidata(app.fh,app);  % store the updated guidata
                %str = sprintf('%g',num);
                %set(cbo,'string',str); % diplay what was stored
                if app.debugging
                    assignin('base','app_ortho',app);
                end
            end
        case 'popupmenu'
            % get the object tag (data struct fieldname)
            tag = get(cbo,'tag');
            tag = tag(1:end-6);  % remove the '_units' ending
            obj = findobj(app.fh,'tag',tag);
            oldunit = get(obj,'Value');
            newunit = get(cbo,'Value');
            if newunit ~= oldunit
                units = get(cbo,'UserData');
                val = app.mat.(tag); % use dynamic fieldname to get value
                val = val * units(newunit);
                str = sprintf('%g',val);
                set(obj,'String',str);
                set(obj,'Value',newunit);
            end
    end
end

function cb_dnum2(cbo,~)
% this callback checks dimensional numeric input and updates app data
    app = guidata(cbo);     % retrieve guidata
    style = get(cbo,'style');  % get the uicontrol style
    switch style
        case 'edit'
            % get the object tag (data struct fieldname)
            tag = get(cbo,'tag');
            % get the units
            obj = findobj(app.fh,'tag',[tag '_units']);
            units = get(obj,'UserData');
            currentunit = get(cbo,'Value');
            % get the stored data associated with the fieldname
            val = app.mat.(tag);  % use dynamic fieldname to get value
            val = val * units(currentunit);  % do unit conversion
            str = get(cbo,'string');  % get the new input
            num = str2double(str);  % try to convert to number
            if isnan(num)
                % input was not valid 
                if isempty(str) % clear value if empty
                    app.mat.(tag) = [];
                    guidata(app.fh,app);  % store the updated guidata
                else % otherwise revert to stored value
                str = sprintf('%g',val);
                set(cbo,'string',str);
                end
            else
                % input valid => save new value
                app.mat.(tag) = num / units(currentunit);
                guidata(app.fh,app);  % store the updated guidata
                %str = sprintf('%g',num);
                %set(cbo,'string',str); % diplay what was stored
                if app.debugging
                    assignin('base','app_ortho',app);
                end
            end
        case 'popupmenu'
            % get the object tag (data struct fieldname)
            tag = get(cbo,'tag');
            tag = tag(1:end-6);  % remove the '_units' ending
            obj = findobj(app.fh,'tag',tag);
            oldunit = get(obj,'Value');
            newunit = get(cbo,'Value');
            if newunit ~= oldunit
                units = get(cbo,'UserData');
                val = app.mat.(tag); % use dynamic fieldname to get value
                val = val * units(newunit);
                str = sprintf('%g',val);
                set(obj,'String',str);
                set(obj,'Value',newunit);
            end
    end
    cb_fccheck(cbo,[]);
end

function cb_propertyset(cbo,~)
    app = guidata(cbo);  % load application data
    selection = get(cbo,'Value');
    set(app.panel,'Visible','off'); % make all invisible
    set(app.panel(selection),'Visible','on') % make selection visible
    if isequal(2,selection)
        set(app.crit,'Visible','on');
    else
        set(app.crit,'Visible','off');
    end
end


function cb_fccheck(cbo,~)
    app = guidata(cbo);  % load application data
    
    fcfields = {'xten','xcmp','yten','ycmp','zten','zcmp','xy','yz','xz',...
        'xycp','yzcp','xzcp','xzit','xzic','yzit','yzic','g1g2',...
        'etal','etat','alp0'};
    fcvalues = cell(numel(fcfields),1);
    
    for kfc = 1:numel(fcfields)
        fcprop = fcfields{kfc};
        fcvalues{kfc} = app.mat.(fcprop);
    end
    
    fcempty = cellfun('isempty',fcvalues);
    
    set(app.crit,'Enable','off');
    if ~any(fcempty(1:9)), set(app.crit(1),'Enable','on'); end  % ESMax
    if ~any(fcempty(1:12)), set(app.crit(2),'Enable','on'); end  % Tsai-Wu
    if ~any(fcempty([1:4,7,13:16])), set(app.crit(3),'Enable','on'); end % Puck
    if ~any(fcempty([1:4,6:7])), set(app.crit(4),'Enable','on'); end % Hashin
    if ~any(fcempty([1:4,7,17:20])), set(app.crit(5),'Enable','on'); end % LaRc
    if ~any(fcempty(1:20)), set(app.crit(6),'Enable','on'); end  % User
    
end


function cb_save(cbo,~,cb_matsave)
    cb_matsave(cbo,[]);  % save the data in the material database
    cb_close(cbo,[]);
end

function cb_close(cbo,~)
    try
        app = guidata(cbo);  % load application data
        position = get(app.fh,'Position');
        settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
        if ~isequal(settings.xy_orthotropic,position(1:2))
            settings.xy_orthotropic = position(1:2);
            writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
        end
        closereq;
    catch ME
        closereq;
        rethrow(ME);
    end
end