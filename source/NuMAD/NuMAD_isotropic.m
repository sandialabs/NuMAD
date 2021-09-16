function NuMAD_isotropic(mat,caller,cb_matsave)
%NUMAD_ISOTROPIC  User interface for isotropic material definition 
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_isotropic launches the interface which allows the user to define
%   the properties of isotropic materials.
%

%===== INITIALIZATION =====================================================
if nargin==3
    if ~isstruct(mat) || ~isfield(mat,'type') || ~strcmp(mat.type,'isotropic')
        errordlg('Argument of this function must be an isotropic material data structure.','Programming Error');
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
    mat.name = 'new isotropic';
    mat.ex = 1e9;
    mat.dens = 1000;
    mat.nuxy = 0.3;
else
    orgname = mat.name;  % remember the original material name in case the user changes it
end

app.debugging = 0;  % with debugging==1, the data structure is sent to the 
                    % workspace after the gui is constructed
%===== GUI CONSTRUCTION ===================================================
units = unitconst();

app.settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
c=cgui([],'init');  c.fig.xy = app.settings.xy_isotropic;
c=cgui(c,'figure','NuMAD - modify isotropic',[500 260],'guiNuMAD-isotropic');
set(c.fig.h,'CloseRequestFcn',@cb_close);

c=cgui(c,'panel','Isotropic material',[10 10 480 240]);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 140     340]);
c=cgui(c,'input','Material Name:',mat.name,'name',@cb_text);
c=cgui(c,'colEdges','',[10 140 250 340]);
c=cgui(c,'dnum','Young''s Modulus (E):',mat.ex,'ex',@cb_dnum,units.stress);
c=cgui(c,'num' ,'Poisson''s Ratio:',mat.nuxy,'nuxy',@cb_num);
c=cgui(c,'dnum','Mass Density:',mat.dens,'dens',@cb_dnum,units.density);
c=cgui(c,'rowHeight','',90);
c=cgui(c,'colEdges','',[10 80 470]);
c=cgui(c,'textbox','Reference:',mat.reference,'reference',@cb_multilinetext);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 340]);
c=cgui(c,'push','Save Changes','','',{@cb_save,cb_matsave});  % note: no @ (already a function handle)
c=cgui(c,'colEdges','',[340 470]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Discard Changes','','',@cb_close);
% jcb: rather that just using 'close', we should maybe check if anything 
% has changed and confirm with user that they wish to discard changes

app.fh = c.fig.h;
app.mat = mat;
app.matname = orgname;  % remember the original material name in case the user changes it
app.caller = caller;
guidata(app.fh,app);  % store the application data as guidata

if app.debugging
    assignin('base','app_iso',app); 
end

end %END GUI CONSTRUCTION (main function)


%===== UTILITY & CALLBACK FUNCTIONS =======================================
function cb_text(cbo,~)
% this callback reads text input and updates app data
    app = guidata(cbo);   % retrieve guidata
    tag = get(cbo,'tag'); % get the object tag (data struct fieldname)
    str = get(cbo,'string');  % get the new input
    app.mat.(tag) = str;
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_iso',app); 
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
        assignin('base','app_iso',app); 
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
            assignin('base','app_iso',app);
        end
    end
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
                    assignin('base','app_iso',app);
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

function cb_save(cbo,~,cb_matsave)
    cb_matsave(cbo,[]);  % save the data in the material database
    cb_close(cbo,[]);
end

function cb_close(cbo,~)
    app = guidata(cbo);  % load application data
    position = get(app.fh,'Position');
    settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
    if ~isequal(settings.xy_isotropic,position(1:2))
        settings.xy_isotropic = position(1:2);
        writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
    end
    closereq;
end
