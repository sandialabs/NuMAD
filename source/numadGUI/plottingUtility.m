function plottingUtility(action)
%PLOTTINGUTILITY  Facilitate plotting of multicolumn data
%   
%   Usage:  plottingUtility()
%     Then enter names of existing workspace variables that contain: 
%        data matrix - with signals in columns and observations in rows
%       signal names - as cell array of strings
%         data units - as cell array of strings
%     Press 'Update' to load X and Y menus
%
%     Use X and Y menus to select columns to plot
%     Use 'Prev' and 'Next' buttons to cycle through Y menu
%     
%     Optionally use the text box to enter other plotting commands.
%     The selected data is available in variables x and y
%     Signal names and units are available in xname, xunit, yname, yunit
%
%   Example:
%     >> A = loadOutData('primary.out')
%     A = 
%            data: [2201x47 double]
%             ffn: [1x86 char]
%         summary: {47x3 cell}
%     >> plottingUtility
%   +---------------------------------+   
%   | Data matrix:   A.data           | 
%   | Signal names:  A.list           | 
%   | Data units:    A.units          | 
%   +---------------------------------+ 
%   
%   See also loadOutData, loadElmData

% ========================================================
%  Written by Jonathan Berg, Sandia National Laboratories
%  v1.0   2009-May-7
%  v1.0.1 2009-May-8   added Data Units to plot labels
%                      added mini command window
%  v1.0.2 2009-May-14  added 'help' documentation 
if ~exist('action','var')
    action = 'initialize';
end

switch action

    case 'initialize'
        gui = findobj('Tag',mfilename);  % Find GUI window
        if isempty(gui)  % If not found, setup new window
            screen = get(0,'ScreenSize'); 
            gui = figure('Color',[0.8 0.8 0.8], ...
                'Position',[0.4*screen(3:4) 400 300], ...
                'IntegerHandle','off', ...
                'MenuBar','none', ...
                'Number','off',...
                'NextPlot','new',...
                'Name','Plotting Utility',...
                'Tag',mfilename);
        else
            figure(gui); clf;  % If found, raise figure and clear
        end
        
        temp = get(gui,'Position');
        fw=temp(3); fh=temp(4);  % get the width and height
        % make label and textbox input for name of data matrix
        uicontrol(...
            'Style','text', ...
            'BackgroundColor',[0.8 0.8 0.8], ...
            'HorizontalAlignment','left', ...
            'Fontsize',10, ...
            'Position',[10 fh-30 80 18], ...
            'String', 'Data matrix: ');
        h.data_matrix = uicontrol(...
            'Style','edit', ...
            'Position',[120 fh-30 200 18], ...
            'String', 'A.data');
        
        % make label and textbox input for signal names variable
        uicontrol(...
            'Style','text', ...
            'BackgroundColor',[0.8 0.8 0.8], ...
            'HorizontalAlignment','left', ...
            'Fontsize',10, ...
            'Position',[10 fh-50 100 18], ...
            'String', 'Signal names: ');
        h.signal_names = uicontrol(...
            'Style','edit', ...
            'Position',[120 fh-50 200 18], ...
            'String', 'A.list');
        
        % make label and textbox input for data units variable
        uicontrol(...
            'Style','text', ...
            'BackgroundColor',[0.8 0.8 0.8], ...
            'HorizontalAlignment','left', ...
            'Fontsize',10, ...
            'Position',[10 fh-70 100 18], ...
            'String', 'Data units: ');
        h.data_units = uicontrol(...
            'Style','edit', ...
            'Position',[120 fh-70 200 18], ...
            'String', 'A.units');
        
        % button which loads signal names into drop-down menus
        uicontrol(...
            'Style','pushbutton', ...
            'Position',[330 fh-70 60 40], ...
            'String', 'Update', ...
            'Callback',[mfilename,'(''update'')']);
        
        % x-variable drop down selector
        uicontrol(...
            'Style','text', ...
            'BackgroundColor',[0.8 0.8 0.8], ...
            'HorizontalAlignment','left', ...
            'Fontsize',10, ...
            'Position',[10 fh-100 20 18], ...
            'String', 'X: ');
        h.xmenu = uicontrol(...
            'Style','popup', ...
            'Position',[30 fh-100 150 20], ...
            'String','no data', ...
            'Fontname','Courier', ...
            'Fontweight','Bold', ...
            'Fontsize',10, ...
            'Value',1, ...
            'Callback',[mfilename,'(''plot'')']);
        
        % y-variable drop down selector
        uicontrol(...
            'Style','text', ...
            'BackgroundColor',[0.8 0.8 0.8], ...
            'HorizontalAlignment','left', ...
            'Fontsize',10, ...
            'Position',[10 fh-125 20 18], ...
            'String', 'Y: ');
        h.ymenu = uicontrol(...
            'Style','popup', ...
            'Position',[30 fh-125 150 20], ...
            'String','no data', ...
            'Fontname','Courier', ...
            'Fontweight','Bold', ...
            'Fontsize',10, ...
            'Value',1, ...
            'Callback',[mfilename,'(''plot'')']);
        
        % Prev and Next buttons
        uicontrol(...
            'Style','pushbutton', ...
            'Position',[30 fh-150 50 20], ...
            'String', '< Prev', ...
            'Callback',[mfilename,'(''previous'')']);
        uicontrol(...
            'Style','pushbutton', ...
            'Position',[130 fh-150 50 20], ...
            'String', 'Next >', ...
            'Callback',[mfilename,'(''next'')']);
        
        % mini comman window
        h.evalmcw = uicontrol(...
            'Style','checkbox',...
            'BackgroundColor',[0.8 0.8 0.8], ...
            'Position',[30 fh-180 70 20], ...
            'String','Evaluate:',...
            'Callback',[mfilename,'(''plot'')']);
        cmdstr = strvcat('[spectra,freq] = pwelch(y);',...
            'Fs = 1 / (x(2)-x(1));',...
            'plot(freq*Fs,10*log10(spectra));',...
            'ylabel([''PSD('' yname '')'']);');
        h.mcw = uicontrol(...
            'Style','edit', ...
            'BackgroundColor',[1 1 1], ...
            'HorizontalAlignment','left', ...
            'Fontname','Courier', ...
            'Position',[30 fh-280 350 100], ...
            'Min',1,'Max',20, ...
            'String', cmdstr);
        
        % save the ui handles
        set(gui,'UserData',h);
        
    case 'update'    
        gui = get(gco,'Parent');  % get handle of GUI from current ui object
        h = get(gui,'UserData');  % recall saved data
        str = get(h.signal_names,'String'); % get user input string
        h.signal_list = evalin('base',str); % evaluate and save result
        str = get(h.data_units,'String');   % get user input string
        h.units_list = evalin('base',str);  % evaluate and save result
        set(h.xmenu,'String',h.signal_list,'Value',1);  % load signal names into the menu, choose item 1
        set(h.ymenu,'String',h.signal_list,'Value',2);  % load signal names into the menu, choose item 2
        
        set(gui,'UserData',h);        % save the signal names for later
        eval([mfilename,'(''plot'')']);  % initialize plot
        
    case 'plot'
        gui = get(gco,'Parent');  % get handle of GUI from current ui object
        h = get(gui,'UserData');  % recall saved data
        datastr = get(h.data_matrix,'String');  % get the Data Matrix string
        xindex = get(h.xmenu,'Value');  % determine x-menu selection
        yindex = get(h.ymenu,'Value');  % determine y-menu selection
        % get the selected columns from the data matrix
        x = evalin('base',sprintf('%s(:,%d)',datastr,xindex));
        y = evalin('base',sprintf('%s(:,%d)',datastr,yindex));
        
        % prepare plot window, create new if necessary
        fig = findobj('Tag','plotwindow');
        if isempty(fig)
            fig = figure('Tag','plotwindow');
        else
            figure(fig); clf;
        end
        xname = strrep(h.signal_list{xindex},'_','\_');
        xunit = regexp(h.units_list{xindex},'[^()]+','match');
        if ~isempty(xunit), xunit = [' (' xunit{:} ')']; end
        yname = strrep(h.signal_list{yindex},'_','\_');
        yunit = regexp(h.units_list{yindex},'[^()]+','match');
        if ~isempty(yunit), yunit = [' (' yunit{:} ')']; end
        if get(h.evalmcw,'Value')
            cmdstr = get(h.mcw,'String');  % load command string
            for k=1:size(cmdstr,1)
                eval(cmdstr(k,:));  % evaluate line by line
            end
        else
            plot(x,y);
            xlabel([xname xunit]);
            ylabel([yname yunit]);
        end
        figure(gui);
        
    case {'next','previous'}
        gui = get(gco,'Parent');  % get handle of GUI from current ui object
        h = get(gui,'UserData');  % recall saved data
        N = numel(h.signal_list);
        yindex = get(h.ymenu,'Value');  % determine y-menu selection
        if (strcmp(action,'next') && yindex<N)
            set(h.ymenu,'Value',yindex+1);
        elseif (strcmp(action,'previous') && yindex>1)
            set(h.ymenu,'Value',yindex-1);
        end
        eval([mfilename,'(''plot'')']);  % update plot
end

