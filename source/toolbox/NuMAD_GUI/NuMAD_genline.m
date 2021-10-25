function NuMAD_genline(cbo,~)
%NUMAD_GENLINE  User interface for blade reference line 
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_genline launches the interface which allows the user to define
%   the blade reference lines for presweep and precurve.
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
units = unitconst();

app.blade = NuMAD_appdata('get','blade');
app.settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
c=cgui([],'init');  c.fig.xy = app.settings.xy_genline;
c=cgui(c,'figure','NuMAD - blade reference line',[800 390],'guiNuMAD-genline');
set(c.fig.h,'CloseRequestFcn',@cb_close);
app.fh = c.fig.h;
%app.bgcolor = c.bgcolor;

c=cgui(c,'panel','Reference Line',[10 10 780 370]);
app.ph = c.panel.h;
c=cgui(c,'rowHeight','',22);
c=cgui(c,'colEdges','',[10 100 240]);
c=cgui(c,'select','Line:',{'Presweep Reference','Precurve Reference'},'selectGenLine',@cb_selectgenline);
c=cgui(c,'select','Method:',{'normal','shear'},'selectMethod',@cb_selectmethod);
c=cgui(c,'select','Line Type:',{'poly','spline','pchip','disabled'},'selectLineType',@cb_selectlinetype);
c=cgui(c,'vpos','next');
c=cgui(c,'colEdges','',[0 10 140]);
c=cgui(c,'select','',{'Insert Row after','Insert Row before','Move down','Move up','Delete'},'rowOperation','');
c=cgui(c,'vpos','previous');
c=cgui(c,'colEdges','',[150 190 240]);
c=cgui(c,'select','Row:',cellstr(num2str(transpose(1:2))),'rowNumber','');
c=cgui(c,'colEdges','',[250 300]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Go','','',@cb_tablerows);


tableColumns = {'Span (m)', 'Offset (m)', 'Slope (m/m)'};
tableData = app.blade.PresweepRef.table;
% for k=1:mat.uniqueLayers
%     layer = mat.layer(k);
%     tableData(k,:) = {layer.layerName layer.thicknessA layer.theta};
% end
app.th=uitable('Parent',c.panel.h,'Position',[10 100 350 120],...
    'ColumnWidth',{100 100 100},'ColumnEditable',[true true true],...
    'Data',tableData,'ColumnName',tableColumns);
% install a callback for data verification
set(app.th,'CellEditCallback',@cb_table);
%set(app.th,'CellSelectionCallback',@cb_rowclick);
rowNumber = get_uictrl(app.fh,'rowNumber');
set(rowNumber.h,'String',cellstr(num2str(transpose(1:size(tableData,1)))),'Value',size(tableData,1));


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
c=cgui(c,'colEdges','',[10 220]);
c=cgui(c,'push','Save Changes','','',@cb_save);
c=cgui(c,'colEdges','',[220 360]);
c=cgui(c,'vpos','previous');
c=cgui(c,'push','Discard Changes','','',@cb_close);

app.ax(1) = axes('Units','normalized','Position',[0.54 0.55 0.44 0.44],'Parent',c.panel.h);
app.ax(2) = axes('Units','normalized','Position',[0.54 0.05 0.44 0.44],'Parent',c.panel.h);

guidata(app.fh,app);  % store the application data as guidata

if app.debugging
    assignin('base','app_comp',app); 
end

obj = findobj(app.fh,'Tag','selectGenLine');
cb_selectgenline(obj,[]);  % initialize the controls

% % display the layer info
% drawlayers(app.fh);
% % select first layer and populate gui controls
% app = guidata(app.fh);
% bdf_selectlayer(app.hlayer(1),[]);

end %END GUI CONSTRUCTION (main function)


%===== UTILITY & CALLBACK FUNCTIONS =======================================
function uictrl = get_uictrl(fh,tag)
    uictrl.h = findobj(fh,'Tag',tag);
    uictrl.string = get(uictrl.h,'String');
    uictrl.value = get(uictrl.h,'Value');
    uictrl.userdata = get(uictrl.h,'UserData');
    if strcmp(get(uictrl.h,'Style'),'popupmenu')
        uictrl.select = uictrl.string{uictrl.value};
    end
end

function update_pp(cbo)
    app = guidata(cbo);
    
    app.blade = calcGenLinePP(app.blade); % update the piecewise polynomial
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_genline',app); 
    end
    
    selectGenLine = get_uictrl(app.fh,'selectGenLine');
    switch selectGenLine.select
        case 'Presweep Reference'
            gl = app.blade.PresweepRef;
        case 'Precurve Reference'
            gl = app.blade.PrecurveRef;
    end
    gl.z = linspace(gl.table(1,1),gl.table(end,1),100);
    gl.offset = ppval(gl.pp,gl.z);
    gl.slope = ppval(gl.dpp,gl.z);   % evaluate derivative spline
    plot(app.ax(1),gl.z,gl.offset,'k-',gl.table(:,1),gl.table(:,2),'ko');
    plot(app.ax(2),gl.z,gl.slope,'k-',gl.table(:,1),gl.table(:,3),'ko');
    ylabel(app.ax(1),'Offset');
    ylabel(app.ax(2),'Slope');
    
end

function cb_selectgenline(cbo,~)
    app = guidata(cbo);
    selectGenLine = get_uictrl(app.fh,'selectGenLine');
    selectMethod = get_uictrl(app.fh,'selectMethod');
    selectLineType = get_uictrl(app.fh,'selectLineType');
    switch selectGenLine.select
        case 'Presweep Reference'
            GenLine = app.blade.PresweepRef;
        case 'Precurve Reference'
            GenLine = app.blade.PrecurveRef;
    end
    set(app.th, 'data', GenLine.table);
    rowNumber = get_uictrl(app.fh,'rowNumber');
    set(rowNumber.h,'String',cellstr(num2str(transpose(1:size(GenLine.table,1)))),'Value',size(GenLine.table,1));
    switch GenLine.method
        case 'normal'
            set(selectMethod.h,'Value',1);
        case 'shear'
            set(selectMethod.h,'Value',2);
    end
    switch GenLine.pptype
        case 'poly'
            set(selectLineType.h,'Value',1);
        case 'spline'
            set(selectLineType.h,'Value',2);
        case 'pchip'
            set(selectLineType.h,'Value',3);
        case 'disabled'
            set(selectLineType.h,'Value',4);
    end 
    
    update_pp(cbo)
end

% function cb_rowclick(cbo,event)
%     app = guidata(cbo);  % load application data
%     obj = findobj(app.fh,'tag','rowNumber');
%     set(obj,'Value',event.Indices(1,1));
% end

function cb_tablerows(cbo,~)
    app = guidata(cbo);  % load application data
    tableData = get(app.th,'data');
    
    obj = findobj(app.fh,'tag','rowOperation');
    rowOperation = get(obj,'Value');
    
    obj = findobj(app.fh,'tag','rowNumber');
    rowNumber = get(obj,'Value');
    
    nRows = size(tableData,1);
    switch rowOperation
        case 1  % insert row after
            tableData = [tableData(1:rowNumber,:);
                         tableData(rowNumber,1) NaN NaN;
                         tableData(rowNumber+1:end,:)]; 
            rowNumberNew = rowNumber+1;
        case 2  % insert row before
            tableData = [tableData(1:rowNumber-1,:);
                         tableData(rowNumber,1) NaN NaN;
                         tableData(rowNumber:end,:)];
            rowNumberNew = rowNumber;
        case 3  % move down
            if rowNumber < nRows
                rowOrder = [1:rowNumber-1 rowNumber+1 rowNumber rowNumber+2:nRows];
                tableData = tableData(rowOrder,:);
                rowNumberNew = rowNumber+1;
            else
                rowNumberNew = rowNumber;
            end
        case 4  % move up
            if rowNumber > 1
                rowOrder = [1:rowNumber-2 rowNumber rowNumber-1 rowNumber+1:nRows];
                tableData = tableData(rowOrder,:);
                rowNumberNew = rowNumber-1;
            else
                rowNumberNew = rowNumber;
            end
        case 5  % delete
            if nRows > 1
                tableData(rowNumber,:) = [];
                rowNumberNew = min(rowNumber,nRows-1);
            else
                rowNumberNew = rowNumber;
            end
    end
    nRows = size(tableData,1);
    set(app.th,'data',tableData);
    obj = findobj(app.fh,'tag','rowNumber');
    set(obj,'String',cellstr(num2str(transpose(1:nRows))));
    set(obj,'Value',rowNumberNew);
    
    selectGenLine = get_uictrl(app.fh,'selectGenLine');
    switch selectGenLine.select
        case 'Presweep Reference'
            app.blade.PresweepRef.table = tableData;
        case 'Precurve Reference'
            app.blade.PrecurveRef.table = tableData;
    end
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_genline',app); 
    end
    
    update_pp(cbo)
end

function cb_selectmethod(cbo,~)
    app = guidata(cbo);
    selectMethod = get_uictrl(app.fh,'selectMethod');
    selectGenLine = get_uictrl(app.fh,'selectGenLine');
    switch selectGenLine.select
        case 'Presweep Reference'
            app.blade.PresweepRef.method = selectMethod.select;
        case 'Precurve Reference'
            app.blade.PrecurveRef.method = selectMethod.select;
    end
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_genline',app); 
    end    
    
    update_pp(cbo)
end

function cb_selectlinetype(cbo,~)
    app = guidata(cbo);
    selectLineType = get_uictrl(app.fh,'selectLineType');
    selectGenLine = get_uictrl(app.fh,'selectGenLine');
    switch selectGenLine.select
        case 'Presweep Reference'
            app.blade.PresweepRef.pptype = selectLineType.select;
        case 'Precurve Reference'
            app.blade.PrecurveRef.pptype = selectLineType.select;
    end
    
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_genline',app); 
    end    
    
    update_pp(cbo)
end

function cb_table(cbo,event)
    app = guidata(cbo);
    tableData = get(cbo, 'data');
    selectGenLine = get_uictrl(app.fh,'selectGenLine');
    switch selectGenLine.select
        case 'Presweep Reference'
            app.blade.PresweepRef.table = tableData;
        case 'Precurve Reference'
            app.blade.PrecurveRef.table = tableData;
    end
    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_genline',app); 
    end    
    
    update_pp(cbo)
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

function cb_save(cbo,~)
    app = guidata(cbo);  % load application data

    blade = NuMAD_appdata('get','blade');
    blade.PresweepRef = app.blade.PresweepRef;
    blade.PrecurveRef = app.blade.PrecurveRef;
    NuMAD_appdata('set','blade',blade);

    guidata(app.fh,app);  % store the updated guidata
    if app.debugging
        assignin('base','app_genline',app); 
    end
    handles = NuMAD_appdata('get','handles');
    main_app = NuMAD_appdata('get','UsedByGUIData_m');
    handles.func.draw_stations(main_app);
    cb_close(cbo,[]);
end

function cb_close(cbo,~)
    app = guidata(cbo);  % load application data
    position = get(app.fh,'Position');
    settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
    if ~isequal(settings.xy_genline,position(1:2))
        settings.xy_genline = position(1:2);
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
