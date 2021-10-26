function AirfoilPrep_v2p2(action)
% AIRFOILPREP_V2P2 GUI implementation of Airfoil Prep v2.2
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%

if nargin==0
    action='initialize';
end

switch action
    case 'initialize'
        OpenGUI;
        
    case 'load_aero'
        [filename,pathname] = uigetfile('*.*');
        if isequal(filename,0)
            % do nothing
        else
            FILEstr=findobj(gcf,'Tag','filename');
            set(FILEstr,'String',filename);
            %don't load old data: ud = get(gcf,'UserData');
            ud.Table2D = LoadAeroData([pathname,filename]);
            ud.fn=filename;
            set(gcf,'UserData',ud);
            
            axes(findobj(gcf,'Tag','Axes1'));
            plot(ud.Table2D(:,1),ud.Table2D(:,2))
            grid on; xlabel('Alpha (deg)');
            legend('CL',2);
            set(gca,'Tag','Axes1');
            
            axes(findobj(gcf,'Tag','Axes2'));
            plot(ud.Table2D(:,1),ud.Table2D(:,3:4))
            grid on; xlabel('Alpha (deg)');
            legend('CD','CM',2);
            set(gca,'Tag','Axes2');
        end
        
    case '3DStall'
        ud = get(gcf,'UserData');
        RPM = str2double(get(findobj(gcf,'Tag','RPM'),'String'));
        R = str2double(get(findobj(gcf,'Tag','R'),'String'));
        V = str2double(get(findobj(gcf,'Tag','V'),'String'));
        rOverR = str2double(get(findobj(gcf,'Tag','rOverR'),'String'));
        Chord = str2double(get(findobj(gcf,'Tag','Chord'),'String'));
        AlphaMinTrend = str2double(get(findobj(gcf,'Tag','Min3D'),'String'));
        AlphaMaxTrend = str2double(get(findobj(gcf,'Tag','Max3D'),'String'));
        AlphaEnd = str2double(get(findobj(gcf,'Tag','AlphaEnd'),'String'));
        [ud.Table3D,ud.Trend3D] = Stall3D(ud.Table2D,RPM,R,V,rOverR,Chord,AlphaMinTrend,AlphaMaxTrend,AlphaEnd);
        set(gcf,'UserData',ud);
        
        axes(findobj(gcf,'Tag','Axes1'));
        plot(ud.Table2D(:,1),ud.Table2D(:,2),'b');
        hold on;
        plot(ud.Table3D(:,1),ud.Table3D(:,2),'x-','Color',[1 0.6 0]);
        plot(ud.Trend3D(:,1),ud.Trend3D(:,2),'+-','LineWidth',2);
        hold off;
        grid on; xlabel('Alpha (deg)'); ylabel('CL');
        b = legend('2-D','3-D','Slope Calc',2);  set(b,'FontSize',8);
        set(gca,'Tag','Axes1');
        
        axes(findobj(gcf,'Tag','Axes2'));
        plot(ud.Table2D(:,1),ud.Table2D(:,3),'b');
        hold on;
        plot(ud.Table3D(:,1),ud.Table3D(:,3),'x-','Color',[1 0.6 0]);
        hold off;
        grid on; xlabel('Alpha (deg)'); ylabel('CD');
        b = legend('2-D','3-D',2);  set(b,'FontSize',8);
        set(gca,'Tag','Axes2');
        
    case 'TableExtrap'
        ud = get(gcf,'UserData');
        CDMax = str2double(get(findobj(gcf,'Tag','CDmax'),'String'));
        UseCM = get(findobj(gcf,'Tag','UseCM'),'Value');
        if get(findobj(gcf,'Tag','FlatPlate'),'Value')
            AlphaFP = str2double(get(findobj(gcf,'Tag','AlphaFP'),'String'));
        else
            AlphaFP = [];
        end
        ud.Table3DE = TableExtrap(ud.Table3D,CDMax,UseCM,AlphaFP);
        set(gcf,'UserData',ud);
        
        axes(findobj(gcf,'Tag','Axes1'));
        plot(ud.Table3DE(:,1),ud.Table3DE(:,2),'b.-',...
            ud.Table3DE(:,1),ud.Table3DE(:,3),'r.-');
        grid on; xlabel('Alpha (deg)');
        b = legend('CL','CD',2);  set(b,'FontSize',8);
        set(gca,'Tag','Axes1');
        
        try
            axes(findobj(gcf,'Tag','Axes2'));
            plot(ud.Table3DE(:,1),ud.Table3DE(:,4),'r.-');
            grid on; xlabel('Alpha (deg)');
            b = legend('CM',1);  set(b,'FontSize',8);
            set(gca,'Tag','Axes2');
        catch
            cla;
        end
        
    case 'DynStall'
        ud = get(gcf,'UserData');
        AlphaMinTrend = str2double(get(findobj(gcf,'Tag','MinStall'),'String'));
        AlphaMaxTrend = str2double(get(findobj(gcf,'Tag','MaxStall'),'String'));
        StallAngle = str2double(get(findobj(gcf,'Tag','StallAngle'),'String'));
        NegStallCn = str2double(get(findobj(gcf,'Tag','NegStallCn'),'String'));
        [ud.StallHeader_ad,ud.StallHeader_wtp,Alpha1,CN,CT,CNSlopeTable,ZeroCN,CNExtrap] = DynStall(ud.Table3D,AlphaMinTrend,AlphaMaxTrend,StallAngle,NegStallCn);
        set(gcf,'UserData',ud);
        
        axes(findobj(gcf,'Tag','Axes2'));
        plot(Alpha1,CN,'m',...
            CNSlopeTable(:,1),CNSlopeTable(:,2),'b.-',...
            [ZeroCN,StallAngle],[0,CNExtrap],'r*-');
        grid on; xlabel('Alpha (deg)');  %xlim([1.5*ZeroCN 1.8*StallAngle]);
        b = legend('CN','Lin. CN','CNStall',2);  set(b,'FontSize',8);
        set(gca,'Tag','Axes2');
        
    case 'print'
        ud = get(gcf,'UserData');
        UseCM = get(findobj(gcf,'Tag','UseCM'),'Value');
        
        % write AeroDyn airfoil file
        fid=fopen(['AD_' ud.fn],'wt');
        for i=1:length(ud.StallHeader_ad)
            fprintf(fid,'%s\n',ud.StallHeader_ad{i,1});
        end
        for k = 1:size(ud.Table3DE,1)
            if (UseCM && size(ud.Table3DE,2)==4)
                fprintf(fid,'%7.2f   %6.3f %7.4f %7.4f\n',ud.Table3DE(k,:));
            else
                fprintf(fid,'%7.2f   %6.3f %7.4f\n',ud.Table3DE(k,:));
            end
        end
        fclose(fid);
        
        % write WT_Perf airfoil file
        fid=fopen(['WTP_' ud.fn],'wt');
        for i=1:length(ud.StallHeader_wtp)
            fprintf(fid,'%s\n',ud.StallHeader_wtp{i,1});
        end
        for k = 1:size(ud.Table3DE,1)
            if (UseCM && size(ud.Table3DE,2)==4)
                fprintf(fid,'%7.2f   %6.3f %7.4f %7.4f\n',ud.Table3DE(k,:));
            else
                fprintf(fid,'%7.2f   %6.3f %7.4f\n',ud.Table3DE(k,:));
            end
        end
        fprintf(fid,'EOT');
        fclose(fid);
end

end %function AirfoilPrep

%==========================================================================
% function that draws user interface
%==========================================================================
function OpenGUI
w = 700;
h = 600;
a = figure('Color',[0.8 0.8 0.8], ...
    'Position',[25 65 w h], ...
    'IntegerHandle','off', ...
    'Number','off',...
    'Name','AeroPrep',...
    'Tag','AeroPrep');

%%%%%%%%%%%%%
% File load %
%%%%%%%%%%%%%
uicontrol('Parent',a,...
    'Style','frame',...
    'String','Input file',...
    'Position',[5 h-95 320 90]);

% control: load button
uicontrol('Parent',a,...
    'Style','pushbutton',...
    'String','Load Airfoil Data',...
    'Position',[40 h-50 100 30],...
    'CallBack',[mfilename,'(''load_aero'')']);

% display loaded filename
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'String','Loaded: ',...
    'Position',[10 h-80 55 17]);

uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'String','',...
    'Position',[55 h-80 200 17],...
    'Tag','filename');

%%%%%%%%%%%
% 3DStall %
%%%%%%%%%%%
h = h - 100;
uicontrol('Parent',a,...
    'Style','frame',...
    'String','3DStall',...
    'Position',[5 h-185 320 180]);

% input: rotor speed
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Rotor speed (rpm)',...
    'Position',[10 h-30 100 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','12',...
    'BackgroundColor',[1 1 1],...
    'Position',[110 h-30 50 17],...
    'Tag','RPM');

% input: rotor radius
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Rotor radius (m)',...
    'Position',[10 h-50 100 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','63',...
    'BackgroundColor',[1 1 1],...
    'Position',[110 h-50 50 17],...
    'Tag','R');

% input: wind speed
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Wind speed (m/s)',...
    'Position',[10 h-70 100 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','12',...
    'BackgroundColor',[1 1 1],...
    'Position',[110 h-70 50 17],...
    'Tag','V');

% input: r/R
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'String','r/R',...
    'Position',[200 h-30 60 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','0.850',...
    'BackgroundColor',[1 1 1],...
    'Position',[260 h-30 50 17],...
    'Tag','rOverR');

% input: chord
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'String','Chord (m)',...
    'Position',[200 h-50 60 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','2.37',...
    'BackgroundColor',[1 1 1],...
    'Position',[260 h-50 50 17],...
    'Tag','Chord');

% input: min/max alpha
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','For CL Slope Calc:',...
    'Position',[200 h-80 120 17]);

uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Min Alpha',...
    'Position',[200 h-100 60 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','-5',...
    'BackgroundColor',[1 1 1],...
    'Position',[260 h-100 50 17],...
    'Tag','Min3D');

uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Max Alpha',...
    'Position',[200 h-120 60 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','5',...
    'BackgroundColor',[1 1 1], ...
    'Position',[260 h-120 50 17],...
    'Tag','Max3D');

% input: alpha end
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','End of full correction:',...
    'Position',[200 h-150 120 17]);

uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Alpha End',...
    'Position',[200 h-170 60 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','30',...
    'BackgroundColor',[1 1 1], ...
    'Position',[260 h-170 50 17],...
    'Tag','AlphaEnd');

% control: 3DStall button
uicontrol('Parent',a,...
    'Style','pushbutton',...
    'String','Calc. 3DStall',...
    'Position',[40 h-150 100 40],...
    'CallBack',[mfilename,'(''3DStall'')']);

%%%%%%%%%%%%%%%
% TableExtrap %
%%%%%%%%%%%%%%%
h = h - 190;
uicontrol('Parent',a,...
    'Style','frame',...
    'String','3DStall',...
    'Position',[5 h-135 320 130]);

% input: aspect ratio
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Aspect ratio',...
    'Position',[10 h-30 70 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','21.7',...
    'BackgroundColor',[1 1 1],...
    'Position',[80 h-30 50 17],...
    'CallBack','set(findobj(gcf,''Tag'',''CDmaxEst''),''String'',num2str(1.11+0.018*str2num(get(gco,''String''))))');
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','right', ...
    'String','1.501',...
    'Position',[140 h-32 50 17],...
    'Tag','CDmaxEst');
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','CDmax estimated from Aspect ratio',...
    'Position',[200 h-50 100 36]);

% input: CDmax
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Max CD',...
    'Position',[10 h-50 70 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','1.46',...
    'BackgroundColor',[1 1 1],...
    'Position',[80 h-50 50 17],...
    'Tag','CDmax');

% checkbox: UseCM
uicontrol('Parent',a,...
    'Style','checkbox',...
    'String','Use CM values',...
    'Value',1,...
    'Position',[20 h-70 100 17],...
    'Tag','UseCM');

% radiobuttions: viterna / flat plate
uicontrol('Parent',a,...
    'Style','radiobutton',...
    'String','Use Viterna method for CD',...
    'Value',1,...
    'Position',[160 h-70 150 17],...
    'Tag','Viterna',...
    'CallBack','set(gco,''Value'',1); set(findobj(gcf,''Tag'',''FlatPlate''),''Value'',0); set(findobj(gcf,''Tag'',''AlphaFP''),''Enable'',''off'')');
uicontrol('Parent',a,...
    'Style','radiobutton',...
    'String','Use flat plate theory for CD',...
    'Value',0,...
    'Position',[160 h-90 150 17],...
    'Tag','FlatPlate',...
    'CallBack','set(gco,''Value'',1); set(findobj(gcf,''Tag'',''Viterna''),''Value'',0); set(findobj(gcf,''Tag'',''AlphaFP''),''Enable'',''on'')');

uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','when abs(angles) greater than               deg',...
    'Value',0,...
    'Position',[176 h-126 130 34]);
uicontrol('Parent',a,...
    'Style','edit',...
    'Enable','off',...
    'String','20',...
    'BackgroundColor',[1 1 1],...
    'Position',[200 h-122 40 17],...
    'Tag','AlphaFP');

% control: extrapolate button
uicontrol('Parent',a,...
    'Style','pushbutton',...
    'String','Extrapolate',...
    'Position',[40 h-120 100 40],...
    'CallBack',[mfilename,'(''TableExtrap'')']);



%%%%%%%%%%%%
% DynStall %
%%%%%%%%%%%%
h = h - 140;
uicontrol('Parent',a,...
    'Style','frame',...
    'String','3DStall',...
    'Position',[5 h-135 320 130]);

% AoA range for CN
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','AoA range for CN slope:',...
    'Position',[10 h-30 120 17]);
% input: min alpha
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Min Alpha',...
    'Position',[10 h-50 70 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','-10',...
    'BackgroundColor',[1 1 1],...
    'Position',[80 h-50 50 17],...
    'Tag','MinStall');
% input: max alpha
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Max Alpha',...
    'Position',[10 h-70 70 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','5',...
    'BackgroundColor',[1 1 1],...
    'Position',[80 h-70 50 17],...
    'Tag','MaxStall');

% input: stall angle
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','Stall angle (deg)',...
    'Position',[170 h-30 80 17]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','14',...
    'BackgroundColor',[1 1 1],...
    'Position',[260 h-30 50 17],...
    'Tag','StallAngle');

% input: CN at stall value for negative AoA
uicontrol('Parent',a,...
    'Style','text',...
    'HorizontalAlignment','left', ...
    'String','CN at stall value for negative AoA',...
    'Position',[170 h-75 100 34]);
uicontrol('Parent',a,...
    'Style','edit',...
    'String','-0.8',...
    'BackgroundColor',[1 1 1],...
    'Position',[260 h-70 50 17],...
    'Tag','NegStallCn');

% control: dynamic stall button
uicontrol('Parent',a,...
    'Style','pushbutton',...
    'String','Calc. DynStall',...
    'Position',[40 h-120 100 40],...
    'CallBack',[mfilename,'(''DynStall'')']);

% control: print button
uicontrol('Parent',a,...
    'Style','pushbutton',...
    'String','Print AD file',...
    'Position',[160 h-120 100 40],...
    'CallBack',[mfilename,'(''print'')']);



%%%%%%%%%
% Plots %
%%%%%%%%%
axes('Parent',a,...
    'Color',[1 1 1],...
    'Position',[.54 .56 .44 .42],...
    'Tag','Axes1');

axes('Parent',a,...
    'Color',[1 1 1],...
    'Position',[.54 .06 .44 .42],...
    'Tag','Axes2');

end % function OpenGUI

%==========================================================================
%==========================================================================
function Table2D = LoadAeroData(fn)
fid = fopen(fn);
if fid==-1
    error('Could not open file: %s',fn)
end
[pn,fn,ext] = fileparts(fn);

switch ext
    case '.forces'  % assume to be ARC2D format
        BadLineMarker = '%';
        
        row = 1;
        Marked = [];
        InData = [];
        while (1)
            line = fgetl(fid);
            if (line == -1)
                break
            elseif (line(1) == BadLineMarker)
                Marked = [Marked row];
                line(1) = [];
            end
            [A, count] = sscanf(line,'%g',[1,inf]);
            if (count ~= size(InData,2) && row>1)
                error('Inconsistent number of data columns')
            end
            InData = [InData; A];
            row = row + 1;
        end
        
        BadLines = false(size(InData,1),1);
        BadLines(Marked) = true;
        Table2D = InData(~BadLines,:);
        %Table2D = InData;
        Table2D(:,2) = [];
        
    case '.dat'  % assume to be XFOIL format
        InData = zeros(0,7);
        for k=1:12
            line = fgetl(fid);
        end
        while (1)
            line = fgetl(fid);
            if (line == -1)
                break
            elseif (isempty(line))
                break
            end
            [A, count] = sscanf(line,'%g',[1,7]);
            InData = [InData; A];
        end
        Table2D = InData(:,[1,2,3,5]);
end

fclose(fid);

end % function LoadAeroData

%==========================================================================
function [Table3D,Trend] = Stall3D(Table2D,RPM,R,V,rOverR,Chord,AlphaMinTrend,AlphaMaxTrend,AlphaEnd,varargin)
% Table3D = stall3d(Table2D,RPM,R,V,rOverR,Chord,AlphaMinTrend,AlphaMaxTrend,AlphaEnd)

% Program to apply Du-Selig and Eggars 3-D stall corrections to 2-D airfoil data
%  C Hansen, Windward Engineering, Dec 2003

a=1; b=1; d=1;  % Basic constants of Selig method.  Do not change except for research
if length(varargin)==1
    abd=varargin{1};
    a=abd(1);
    b=abd(2);
    d=abd(3);
end

DToR = pi/180;

% Call Update2DSlope
nTable2D = size(Table2D,1);
Alpha2D = Table2D(:,1);
CL2D = Table2D(:,2);
CD2D = Table2D(:,3);
rows = ((Alpha2D >= AlphaMinTrend) & (Alpha2D <= AlphaMaxTrend));
TrendInput = Table2D(rows,1:2);
P = polyfit(TrendInput(:,1),TrendInput(:,2),1);
CLSlope = P(1);
CLIntercept = P(2);
AlphaZero = -CLIntercept/CLSlope;
Trend = [TrendInput(:,1),polyval(P,TrendInput(:,1))];

% calculate and display some constants
% uses original notation for Selig variables
coverr = Chord / R / rOverR;  % c/r
Lambda = RPM * pi/30 * R / sqrt(V^2 + (RPM * pi/30 * R)^2);
Expon = d / Lambda / rOverR;
FL = 1 / (CLSlope / DToR) * ((1.6 * coverr / 0.1267)*(a - coverr ^ Expon) / (b + coverr ^ Expon) - 1);

% calculate 3D values
AlphaRad = DToR * Alpha2D;
CN2D = CL2D .* cos(AlphaRad) + CD2D .* sin(AlphaRad);
CT2D = CL2D .* sin(AlphaRad) - CD2D .* cos(AlphaRad);

CLP = CLSlope * (Alpha2D - AlphaZero);
DelCL1 = FL * (CLP - CL2D);

adj = ones(nTable2D,1);
rows = (Alpha2D > AlphaEnd);
adj(rows) = ((90 - Alpha2D(rows)) ./ (90 - AlphaEnd)).^2;

CL3D = CL2D + DelCL1 .* adj;
DelCL2 = CL3D - CL2D;

DelCD = DelCL2 .* (sin(AlphaRad) - 0.12 * cos(AlphaRad)) ./ (cos(AlphaRad) + 0.12 * sin(AlphaRad));
CD3D = CD2D + DelCD;

% Display the results if in the appropriate range
%  (note this method does not work for very large or small angles)
Table3D = [Alpha2D, CL3D, CD3D, Table2D(:,4)];
rows = ((Alpha2D >= -10) & (Alpha2D <= 90));
Table3D = Table3D(rows,:);

end


%==========================================================================
function OutputTable = TableExtrap(InputTable,CDMax,UseCM,varargin)
% OutputTable = TableExtrap(InputTable,CDMax,UseCM,[AlphaFP])

%Foilcheck program converted to Visual Basic routine
% CM not included in this version!
% CH Windward Engineering, Dec 2003

UseFlatPlate = false;
if length(varargin)==1
    %Get angle when flap plate theory will be applied (if that option selected)
    AlphaFP = varargin{1};
    UseFlatPlate = true;
end

DToR = pi/180;

% read in original table values
nTable1 = size(InputTable,1);
Alpha1 = InputTable(:,1);
CL1 = InputTable(:,2);
CD1 = InputTable(:,3);
if (UseCM)
    CM1 = InputTable(:,4);
else
    CM1 = zeros(size(Alpha1));
end

% apply viterna extrapolation
% ==============================
VAlphaLo = Alpha1(1);
VAlphaHi = Alpha1(end);

% Find CL and CD at VAlphaHi, the upper matching point
for i = 1:nTable1
    if (abs(Alpha1(i) - VAlphaHi) < 0.01)
        CLRef = CL1(i);
        CDRef = CD1(i);
        %fprintf(' CLRef = %g\n',CLRef);
        break
    end
    if (i == nTable1)  % we shouldn't get here
        error('Error: Your upper matching point for the Viterna calculation was not found in the original airfoil table');
    end
end

% get Viterna coefficients
SAlpha = sin(VAlphaHi * DToR);
CAlpha = cos(VAlphaHi * DToR);

A2 = (CLRef - CDMax / 2 * sin(2 * VAlphaHi * DToR)) * SAlpha / CAlpha ^ 2;
B2 = (CDRef - CDMax * SAlpha ^ 2) / CAlpha;

% get pitching moment coefficients (if needed)
if (UseCM)
    % call CMCoeff
    hi = find(Alpha1 >= VAlphaHi,1);
    CLHi = CL1(hi);
    CDHi = CD1(hi);
    CMHi = CM1(hi);
    
    FoundZeroLift = false;
    % get CM at angle of zero lift (CM0)
    for i = 1:nTable1-1
        if (abs(Alpha1(i)) < 20 && CL1(i) <= 0 && CL1(i+1) >= 0)
            p = -CL1(i) / (CL1(i+1) - CL1(i));
            CM0 = CM1(i) + p*(CM1(i+1) - CM1(i));
            FoundZeroLift = true;
            break
        end
    end
    if (~FoundZeroLift)  % zero lift not in range of orig table, use first two points
        p = -CL1(1) / (CL1(2) - CL1(1));
        CM0 = CM1(1) + p*(CM1(2) - CM1(1));
    end
    
    XM = (-CMHi + CM0) / (CLHi * cos(VAlphaHi * DToR) + CDHi * sin(VAlphaHi * DToR));
    CMCoef = (XM - 0.25) / tan((VAlphaHi - 90) * DToR);
    
end

% call NewTable
Step = 10;   % angle of attack step in deg
SizeMax = 200;  % dimension of Alpha arrays
HiEnd = 0;
remainder = VAlphaHi - floor(VAlphaHi / Step) * Step;

Alpha2 = zeros(SizeMax,1);
CL2 = zeros(SizeMax,1);
CD2 = zeros(SizeMax,1);
CM2 = zeros(SizeMax,1);

Alpha2(1) = -180;
count = 1;
for i = 2:SizeMax
    Alpha2(i) = Alpha2(i-1) + Step;
    if (Alpha2(i) > 180 - Step)         % stop when we reach 180 deg
        Alpha2(i) = 180;
        nTable2 = i;
        break
    end
    if (Alpha2(i) >= VAlphaLo)
        if (Alpha2(i) - Step < VAlphaHi) %in the range of the original table, use those alphas
            Alpha2(i) = Alpha1(count);    %save these values and write to worksheet
            CL2(i) = CL1(count);
            CD2(i) = CD1(count);
            CM2(i) = CM1(count);
            count = count + 1;    % number of original alpha values we've used
        elseif (HiEnd == 0)   % we just left the upper end of the original alpha table
            HiEnd = 1;
            Alpha2(i) = VAlphaHi - remainder + Step; %find first multiple of 10 after ValphaHi
        end
    end
end

% call ViternaFill
CLAdj = 0.7;  % adjustment factor for CL when abs(Alpha)>90

lo = find(Alpha1 >= VAlphaLo,1);
hi = find(Alpha1 >= VAlphaHi,1);
CLLo = InputTable(lo,2);
CDLo = InputTable(lo,3);
CLHi = InputTable(hi,2);
CDHi = InputTable(hi,3);

%Get angle when flap plate theory will be applied (if that option selected)
%%AlphaFP = Cells(11, 2).Value

for i = 1:nTable2    %step through all alphas.  Notice that if angle is
    %within original table range the value is not changed
    Alfa = Alpha2(i);
    SAlfa = sin(Alfa * DToR);
    CAlfa = cos(Alfa * DToR);
    if (Alfa > 180 || Alfa < -180)
        error(['Angle of attack = ',Alfa,' outside range + to -180 deg in Viterna calculation'])
    elseif (Alfa >= VAlphaHi && Alfa <= 90)
        CL2(i) = CDMax / 2 * sin(2 * Alfa * DToR) + A2 * CAlfa ^ 2 / SAlfa;
        CD2(i) = CDMax * SAlfa ^ 2 + B2 * CAlfa;
    elseif (Alfa > 90 && Alfa <= 180 - VAlphaHi)
        Ang = 180 - Alfa;
        Sang = sin(Ang * DToR);
        Cang = cos(Ang * DToR);
        CL2(i) = CLAdj * (-CDMax / 2 * sin(2 * Ang * DToR) - A2 * Cang ^ 2 / Sang);
        CD2(i) = CDMax * Sang ^ 2 + B2 * Cang;
    elseif (Alfa > 180 - VAlphaHi && Alfa <= 180)
        Ang = Alfa - 180;
        Sang = sin(Ang * DToR);
        Cang = cos(Ang * DToR);
        CL2(i) = Ang / VAlphaHi * CLHi * CLAdj;
        CD2(i) = CDMax * Sang ^ 2 + B2 * Cang;
    elseif (Alfa > -VAlphaHi && Alfa < VAlphaLo)
        CL2(i) = CLAdj * (-CLHi + (Alfa + VAlphaHi) / (VAlphaHi + VAlphaLo) * (CLHi + CLLo));
        CD2(i) = CDLo + (-Alfa + VAlphaLo) * (CDHi - CDLo) / (VAlphaHi + VAlphaLo);
    elseif (Alfa <= -VAlphaHi && Alfa >= -90)
        Ang = -Alfa;
        Sang = sin(Ang * DToR);
        Cang = cos(Ang * DToR);
        CL2(i) = CLAdj * (-CDMax / 2 * sin(2 * Ang * DToR) - A2 * Cang ^ 2 / Sang);
        CD2(i) = CDMax * Sang ^ 2 + B2 * Cang;
    elseif (Alfa < -90 && Alfa >= -180 + VAlphaHi)
        Ang = 180 + Alfa;
        Sang = sin(Ang * DToR);
        Cang = cos(Ang * DToR);
        CL2(i) = CLAdj * (CDMax / 2 * sin(2 * Ang * DToR) + A2 * Cang ^ 2 / Sang);
        CD2(i) = CDMax * Sang ^ 2 + B2 * Cang;
    elseif (Alfa < -180 + VAlphaHi && Alfa >= -180)
        Ang = 180 + Alfa;
        Sang = sin(Ang * DToR);
        Cang = cos(Ang * DToR);
        CL2(i) = Ang / VAlphaHi * CLHi * CLAdj;
        CD2(i) = CDMax * Sang ^ 2 + B2 * Cang;
    end % if
    
    if (CD2(i) < 0)   %watch out for negative CD's
        CD2(i) = 0.01;
    end
end

%Change to using CD from flat plate theory if desired
if (UseFlatPlate)
    for i = 1:nTable2
        Alfa = Alpha2(i);
        if (abs(Alfa) >= AlphaFP)
            CD2(i) = CL2(i) * tan(Alfa * DToR);
        end
    end
end

%Do CM calculations and write value to output table
if (UseCM)
    for i = 1:nTable2
        Alfa = Alpha2(i);
        if (Alfa >= VAlphaLo && Alfa <= VAlphaHi)
            continue %no action needed
        end
        if (Alfa > -165 && Alfa < 165)
            if (abs(Alfa) < 0.01)
                CM2(i) = CM0;
            else
                if (Alfa > 0)
                    x = CMCoef * tan((Alfa - 90) * DToR) + 0.25; % coefficient used in CM calc
                    CM2(i) =   CM0 - x * (CL2(i) * cos(Alfa * DToR) + CD2(i) * sin(Alfa * DToR));
                else
                    % correction by Scott Larwood for negative Alfa (AirfoilPrep_v2.2.1)
                    x = CMCoef * tan((-Alfa - 90) * DToR) + 0.25; % coefficient used in CM calc
                    CM2(i) = -(CM0 - x * (-CL2(i) * cos(-Alfa * DToR) + CD2(i) * sin(-Alfa * DToR)));
                end
            end
        else
            switch (Alfa)
                case 165
                    CM2(i) = -0.4;
                case 170
                    CM2(i) = -0.5;
                case 175
                    CM2(i) = -0.25;
                case 180
                    CM2(i) = 0;
                case -165
                    CM2(i) = 0.35;
                case -170
                    CM2(i) = 0.4;
                case -175
                    CM2(i) = 0.2;
                case -180
                    CM2(i) = 0;
                otherwise
                    error('Angle encountered for which there is no CM table value (near +/-180°).  Program will stop');
            end
        end
    end %for i
    
    OutputTable = [Alpha2(1:nTable2), CL2(1:nTable2), CD2(1:nTable2), CM2(1:nTable2)];
    
else
    
    OutputTable = [Alpha2(1:nTable2), CL2(1:nTable2), CD2(1:nTable2)];
    
end %if (UseCM)

end % function


%==========================================================================
function [StallHeader_ad,StallHeader_wtp,Alpha1,CN,CT,CNSlopeTable,ZeroCN,CNExtrap] = DynStall(InputTable,AlphaMinTrend,AlphaMaxTrend,StallAngle,NegStallCn)
% StallHeader = DynStall(InputTable,AlphaMinTrend,AlphaMaxTrend,StallAngle,NegStallCn)

%Foilcheck program converted to Visual Basic routine
% CM not included in this version!
% CH Windward Engineering, Dec 2003

DToR = pi/180;

Alpha1 = InputTable(:,1);
CL1 = InputTable(:,2);
CD1 = InputTable(:,3);
CM1 = InputTable(:,4);
CN = CL1 .* cos(DToR * Alpha1) + CD1 .* sin(DToR * Alpha1);
CT = CL1 .* sin(DToR * Alpha1) - CD1 .* cos(DToR * Alpha1);

% get CN slope
rows = ((Alpha1 >= AlphaMinTrend) & (Alpha1 <= AlphaMaxTrend));
CNSlopeTable = [Alpha1(rows),CN(rows)];
P = polyfit(CNSlopeTable(:,1),CNSlopeTable(:,2),1);
CNSlope = P(1)*180/pi;           % Cn slope for zero lift
ZeroCN = -P(2)/CNSlope*180/pi;   % Zero Cn angle of attack (deg)
CNExtrap = CNSlope*pi/180*(StallAngle-ZeroCN);

% get angle of attack for CDmin (search only in range -20 to +20 deg)
rows = find(((Alpha1 > -20) & (Alpha1 < 20)));
[y,i] = min(CD1(rows));
CDMin = CD1(rows(i));
AlphaCDMin = Alpha1(rows(i));

StallHeader_ad{1,1} = sprintf('title line 1');
StallHeader_ad{2,1} = sprintf('title line 2');
StallHeader_ad{3,1} = sprintf('%4i       Number of airfoil tables in this file',1);
StallHeader_ad{4,1} = sprintf('%7.2f    Table ID parameter',0);
StallHeader_ad{5,1} = sprintf('%7.2f    Stall angle (deg)',StallAngle);
StallHeader_ad{6,1} = sprintf('%7.2f    No longer used, enter zero',0.0);
StallHeader_ad{7,1} = sprintf('%7.2f    No longer used, enter zero',0.0);
StallHeader_ad{8,1} = sprintf('%7.2f    No longer used, enter zero',0.0);
StallHeader_ad{9,1} = sprintf('%9.4f  Zero Cn angle of attack (deg)',ZeroCN);
StallHeader_ad{10,1} = sprintf('%9.4f  Cn slope for zero lift (dimensionless)',CNSlope);
StallHeader_ad{11,1} = sprintf('%9.4f  Cn extrapolated to value at positive stall angle of attack',CNExtrap);
StallHeader_ad{12,1} = sprintf('%9.4f  Cn at stall value for negative angle of attack',NegStallCn);
StallHeader_ad{13,1} = sprintf('%7.2f    Angle of attack for minimum CD (deg)',AlphaCDMin);
StallHeader_ad{14,1} = sprintf('%9.4f  Minimum CD value',CDMin);

StallHeader_wtp{1,1} = sprintf('AeroDyn airfoil file.  Compatible with AeroDyn v13.0.');
StallHeader_wtp{2,1} = sprintf('title line 2');
StallHeader_wtp{3,1} = sprintf('title line 3');
StallHeader_wtp{4,1} = sprintf('%4i       Number of airfoil tables in this file',1);
StallHeader_wtp{5,1} = sprintf('%7.2f    Table ID parameter',0);
StallHeader_wtp{6,1} = sprintf('%7.2f    Stall angle (deg)',StallAngle);
StallHeader_wtp{7,1} = sprintf('%9.4f  Zero Cn angle of attack (deg)',ZeroCN);
StallHeader_wtp{8,1} = sprintf('%9.4f  Cn slope for zero lift (dimensionless)',CNSlope);
StallHeader_wtp{9,1} = sprintf('%9.4f  Cn extrapolated to value at positive stall angle of attack',CNExtrap);
StallHeader_wtp{10,1} = sprintf('%9.4f  Cn at stall value for negative angle of attack',NegStallCn);
StallHeader_wtp{11,1} = sprintf('%7.2f    Angle of attack for minimum CD (deg)',AlphaCDMin);
StallHeader_wtp{12,1} = sprintf('%9.4f  Minimum CD value',CDMin);
end
