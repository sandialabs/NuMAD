function NuMAD_composite(mat,caller,cb_matsave)
%NUMAD_COMPOSITE  NuMAD composite material user interface
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_composite launches the interface which allows the user to
%   define the layers of a composite material.



%===== CREDITS & CHANGELOG ================================================
%Developed by Wind & Water Power Technologies, Sandia National Laboratories
%2011.01.04  JCB: first draft
%2011.01.28  JCB: fully functional data input using uitable
%     ToDo - visualization of composite layers

if nargin==3
    if ~isstruct(mat) || ~isfield(mat,'type') || ~strcmp(mat.type,'composite')
        errordlg('Argument of this function must be a composite material data structure.','Programming Error');
        return;
    end
else
    errordlg('Script must be called from main materials window.','Programming Error');
    return;
end

callapp = guidata(caller);
app.userpath = callapp.userpath;
isoortho = callapp.isotropic | callapp.orthotropic;
app.isoortho = callapp.matlist(isoortho);

if numel(app.isoortho) < 1
    error('At least one Isotropic or Orthotropic/Layer material must be defined to create a composite material.');
end

if strcmp(mat.name,'__NewMatFlag')
    orgname = mat.name;  % remember the original material name in case the user changes it
    mat.name = 'new composite';
    mat.reference = '';
    mat.thicknessType = 'Constant';
    mat.uniqueLayers = 1;
    mat.symmetryType = 'none';
    mat.layer(1).layerName = app.isoortho{1};
    mat.layer(1).thicknessA = 0.001;
    mat.layer(1).thicknessB = 0.001;
    mat.layer(1).quantity = 1;
    mat.layer(1).theta = 0;
else
    orgname = mat.name;  % remember the original material name in case the user changes it
end    


app.debugging = 0;  % with debugging==1, the data structure is sent to the 
                    % workspace after the gui is constructed
%===== GUI CONSTRUCTION ===================================================
units = unitconst();

app.settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
c=cgui([],'init');  c.fig.xy = app.settings.xy_composite;
c=cgui(c,'figure','NuMAD - modify composite',[550 490],'guiNuMAD-composite');
set(c.fig.h,'CloseRequestFcn',@cb_close);
app.fh = c.fig.h;
%app.bgcolor = c.bgcolor;

c=cgui(c,'panel','Composite material',[10 10 530 470]);
app.ph = c.panel.h;
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 100 515]);
c=cgui(c,'input','Material Name:',mat.name,'name',@cb_text);
c=cgui(c,'rowHeight','',90);
c=cgui(c,'colEdges','',[10 100 515]);
c=cgui(c,'textbox','Reference:',mat.reference,'reference',@cb_multilinetext);
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 100 160]);
c=cgui(c,'select','Symmetry:',{'none','even','odd'},'symmetryType',@cb_symmetryType);
initialize(app.fh,'symmetryType',mat.symmetryType);
c=cgui(c,'vpos','next');
c=cgui(c,'colEdges','',[0 10 140]);
c=cgui(c,'select','',{'Insert row after','Insert row before','Move down','Move up','Delete'},'rowOperation','');
c=cgui(c,'vpos','previous');
c=cgui(c,'colEdges','',[150 190 240]);
c=cgui(c,'select','Row:',cellstr(num2str(transpose(1:mat.uniqueLayers))),'rowNumber','');
c=cgui(c,'colEdges','',[250 300]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Go','','',@cb_editrows);

tableColumns = {'Material', 'Thickness|(mm)', 'Qty.', 'Orientation|(degrees)'};
tableData = cell(mat.uniqueLayers,4);
for k=1:mat.uniqueLayers
    layer = mat.layer(k);
    % multiply thickness by 1000 to convert to mm
    tableData(k,:) = {layer.layerName layer.thicknessA*1e3 layer.quantity layer.theta};
end
app.th=uitable('Parent',c.panel.h,'Position',[10 40 500 220],...
    'ColumnWidth',{240 80 40 80},'ColumnEditable',[true true true true],...
    'Data',tableData,'ColumnName',tableColumns);
% restrict set of choices for material column
set(app.th,'ColumnFormat',{app.isoortho [] [] []}); 
% install a callback for data verification
set(app.th,'CellEditCallback',@cb_layertable);

% c=cgui(c,'colEdges','',[10 50 80]);
% c=cgui(c,'num','Layer:',1,'layer','');
% set(c.ctrl.h,'Enable','off');
% c=cgui(c,'vpos','previous');
% c=cgui(c,'colEdges','',[100 150 250]);
% c=cgui(c,'select','Material:',app.isoortho,'layerName','');
% c=cgui(c,'vpos','previous');
% c=cgui(c,'colEdges','',[270 330 420]);
% c=cgui(c,'input','Thickness:','','thickness','');
% c=cgui(c,'vpos','previous');
% c=cgui(c,'colEdges','',[440 500 570]);
% c=cgui(c,'input','Orientation:','','theta','');

c=cgui(c,'vpos','value',35);
c=cgui(c,'colEdges','',[10 340]);
c=cgui(c,'push','Save Changes','','',{@cb_save,cb_matsave});  % note: no @ (already a function handle)
c=cgui(c,'colEdges','',[340 470]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Discard Changes','','',@cb_close);

app.mat = mat;
app.matname = orgname;  % remember the original material name in case the user changes it
app.caller = caller;
guidata(app.fh,app);  % store the application data as guidata

if app.debugging
    assignin('base','app_comp',app); 
end

% % display the layer info
% drawlayers(app.fh);
% % select first layer and populate gui controls
% app = guidata(app.fh);
% bdf_selectlayer(app.hlayer(1),[]);

end %END GUI CONSTRUCTION (main function)


%===== UTILITY & CALLBACK FUNCTIONS =======================================
function cb_text(cbo,~)
% this callback reads text input and updates app data
    app = guidata(cbo);  % load application data
    tag = get(cbo,'tag'); % get the object tag (data struct fieldname)
    str = get(cbo,'string');  % get the new input
    app.mat.(tag) = str;
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_comp',app); 
    end
end

function cb_multilinetext(cbo,~)
% this callback reads text input and updates app data
    app = guidata(cbo);  % load application data
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
        assignin('base','app_comp',app); 
    end
end

function cb_symmetryType(cbo,~)
% this callback processes changes to the symmetry dropdown menu
    app = guidata(cbo);   % retrieve application data
    menulist = get(cbo,'string');
    selectval = get(cbo,'value');
    selectstr = menulist{selectval};
    app.mat.symmetryType = selectstr;
    
    guidata(app.fh,app);  % store the updated application data
    if app.debugging
        assignin('base','app_comp',app);
    end
end

function cb_editrows(cbo,~)
    app = guidata(cbo);  % load application data
    tableData = get(app.th,'data');
    
    obj = findobj(app.fh,'tag','rowOperation');
    rowOperation = get(obj,'Value');
    
    obj = findobj(app.fh,'tag','rowNumber');
    rowNumber = get(obj,'Value');
    
    nLayers = size(tableData,1);
    switch rowOperation
        case 1  % insert layer after
            newLayer = tableData(rowNumber,:);
            tableData = [tableData(1:rowNumber,:);
                         newLayer;
                         tableData(rowNumber+1:end,:)]; 
            rowNumberNew = rowNumber+1;
        case 2  % insert layer before
            newLayer = tableData(rowNumber,:);
            tableData = [tableData(1:rowNumber-1,:);
                         newLayer;
                         tableData(rowNumber:end,:)];
            rowNumberNew = rowNumber;
        case 3  % move down
            if rowNumber < nLayers
                layerOrder = [1:rowNumber-1 rowNumber+1 rowNumber rowNumber+2:nLayers];
                tableData = tableData(layerOrder,:);
                rowNumberNew = rowNumber+1;
            else
                rowNumberNew = rowNumber;
            end
        case 4  % move up
            if rowNumber > 1
                layerOrder = [1:rowNumber-2 rowNumber rowNumber-1 rowNumber+1:nLayers];
                tableData = tableData(layerOrder,:);
                rowNumberNew = rowNumber-1;
            else
                rowNumberNew = rowNumber;
            end
        case 5  % delete
            if nLayers > 1
                tableData(rowNumber,:) = [];
                rowNumberNew = min(rowNumber,nLayers-1);
            else
                warndlg('There must be at least one layer in a composite material.','Operation Not Permitted');
                rowNumberNew = rowNumber;
            end
    end
    nLayers = size(tableData,1);
    set(app.th,'data',tableData);
    obj = findobj(app.fh,'tag','rowNumber');
    set(obj,'String',cellstr(num2str(transpose(1:nLayers))));
    set(obj,'Value',rowNumberNew);
    
%     guidata(app.fh,app);  % store the updated guidata
%     if app.debugging
%         assignin('base','app_comp',app); 
%     end
end

function cb_layertable(cbo,event)
%     app = guidata(cbo);  % load application data
%     tableData = get(cbo, 'data');
%     app.tableData = tableData;
    
%     guidata(app.fh,app);  % store the updated guidata
%     if app.debugging
%         assignin('base','app_comp',app); 
%     end
end

% jcb: delete the following function block later
% This block shows how to write a callback for a uitable
% function AgeVerificationCallback(o, e)
%    if (e.Indices(2) == 2 && ...
%        (e.NewData < 0 || e.NewData > 120))
%        tableData = get(o, 'data');
%        tableData{e.Indices(1), e.Indices(2)} = e.PreviousData;
%        set(o, 'data', tableData);
%        error('Age value must be between 0 and 120.')
%    end
% end

function cb_save(cbo,~,cb_matsave)
    app = guidata(cbo);  % load application data
    tableData = get(app.th, 'data');

    nRows = size(tableData,1);
    app.mat.uniqueLayers = nRows;
    for k=1:nRows
        app.mat.layer(k).layerName = tableData{k,1};
        % divide thickness by 1000 to convert back to m
        app.mat.layer(k).thicknessA = tableData{k,2}*1e-3;
        app.mat.layer(k).thicknessB = tableData{k,2}*1e-3;
        app.mat.layer(k).quantity = tableData{k,3};
        app.mat.layer(k).theta = tableData{k,4};
    end

    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_comp',app); 
    end
    cb_matsave(cbo,[]);  % save the data in the material database
    cb_close(cbo,[]);
end

function cb_close(cbo,~)
    app = guidata(cbo);  % load application data
    position = get(app.fh,'Position');
    settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
    if ~isequal(settings.xy_composite,position(1:2))
        settings.xy_composite = position(1:2);
        writeNuMADsettings(settings,fullfile(app.userpath,'settings.txt'));
    end
    closereq;
end


function initialize(fh,tag,value)
    obj = findobj(fh,'tag',tag);
    style = get(obj,'style');
    switch style
        case 'popupmenu'
            menulist = get(obj,'string');
            k=find(strcmp(value,menulist)==1);
            set(obj,'value',k);
        otherwise
            errordlg('The initialize() function does not currently handle that type of uicontrol.','Programming Error');
            return;
    end
end

% function drawlayers(cbo,~)
%     app = guidata(cbo);
%     mat = app.mat;
%     
%     % old NuMAD:  red, blue, green, brown, yellow, light blue, orange,
%     %             dark green, violet, light brown, cyan, red, light green,
%     %             dark red, light orange, ...
%     
%     %lines(7);
%     layercolors = [0.00 0.00 1.00;
%                    0.00 0.50 0.00;
%                    1.00 0.00 0.00;
%                    0.00 0.75 0.75;
%                    0.75 0.00 0.75;
%                    0.75 0.75 0.00;
%                    0.25 0.25 0.25];
%     Ncolors = size(layercolors,1);
%     
%     if isfield(app,'hlayer')
%         delete(app.hlayer);
%     end
%     strpos = -10;
%     for k = 1:numel(mat.layer)
%        str = sprintf('Layer %d, %s, %g, %g',...
%            k, mat.layer(k).layerName, mat.layer(k).thicknessA, mat.layer(k).theta);
%        strpos = strpos + 20;
%        lyrclr = layercolors(mod(k-1,Ncolors)+1,:);  % loop thru colors
%        h(k) = uicontrol(app.fh,...
%            'Style','text',...
%            'String',str,...
%            'Value',k,...
%            'FontSize',10,...
%            'HorizontalAlignment','left',...
%            'ForegroundColor',lyrclr,...
%            'Position',[10 strpos 200 20],...
%            'ButtonDownFcn',@bdf_selectlayer);
%     end
%     app.hlayer = h;
%     % store application data
%     guidata(app.fh,app);  
%     if app.debugging
%       assignin('base','app_comp',app);
%     end
% end
% 
% function bdf_selectlayer(cbo,~)
% % this "button down function" responds to mouse clicks on graphics
% app = guidata(cbo);  % retrieve the application data
% %click = get(gcbf,'SelectionType');  % find out which mouse button was pressed
% %switch click
% %    case {'normal','extend','alt'}  % any mouse button
%         % look for active (already selected) layer
%         %active = findobj(gcbf,'Tag','active layer');
%         active = findobj(gcf,'Tag','active layer');
%         if ~isempty(active)
%             % layer found, deselect the layer
%             set(active,'FontWeight','normal','Tag','',...
%                 'BackgroundColor',app.bgcolor);
%         end
%         % make the cbo the active layer
%         set(cbo,'FontWeight','bold','Tag','active layer',...
%             'BackgroundColor',[1 1 1]);
%         readActiveLayer(cbo);
%         
% %    case 'open'  % double mouse click (any button)
% %end
% end
% 
% function readActiveLayer(cbo)
% % this function populates the gui controls with the active layer's
% % parameter values
%     % load application data
%     app = guidata(cbo);
%     % determine which layer is active
%     k = get(cbo,'Value');
%     
%     tag = 'layer';
%     obj = findobj(app.fh,'Tag',tag);
%     val = k;
%     str = sprintf('%g',val);
%     set(obj,'string',str);
% 
%     tag = 'layerName';
%     obj = findobj(app.fh,'Tag',tag);
%     list = get(obj,'String');
%     n = strcmp(app.mat.layer(k).(tag),list);
%     set(obj,'Value',find(n==1));
%     
%     tag = 'thickness';
%     obj = findobj(app.fh,'Tag',tag);
%     val = app.mat.layer(k).('thicknessA');
%     str = sprintf('%g',val);
%     set(obj,'string',str);
%     
%     tag = 'theta';
%     obj = findobj(app.fh,'Tag',tag);
%     val = app.mat.layer(k).(tag);
%     str = sprintf('%g',val);
%     set(obj,'string',str);
%     
% end
