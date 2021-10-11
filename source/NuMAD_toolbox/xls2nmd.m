function varargout = xls2nmd(varargin)
%XLS2NMD  Convert MS Excel file NuMAD.xls into NuMAD input files
% **********************************************************************
% *           Part of the SNL Wind Turbine Analysis Toolbox            *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% **********************************************************************
%   xls2nmd(noinput)
%   Reads NuMAD.xls (or *.xlsx) Excel file and creates NMD input file and
%   MatDBsi.txt material input file.  NMD file name is determined by entry
%   in cell 1,1 of the "Geometry" sheet in the Excel worksheet.
%

%===== CREDITS & CHANGELOG ================================================
% 2011.06.06 brr: creation
% yyyy.mm.dd initials: description

if nargin>0
    xlsname = varargin{1};
    if ~exist(xlsname,'file')
        error('File "%s" not found.',xlsname);
    end
else
    if exist('NuMAD.xlsx','file')
        xlsname = 'NuMAD.xlsx';
    elseif exist('NuMAD.xls','file')
        xlsname = 'NuMAD.xls';
    else
        disp('"NuMAD.xls(x)" not found. Please specify .xls file.')
        return
    end
end

% Read in data from Excel file, NuMAD.xlsx
try
    [~,sheets]=xlsfinfo(xlsname);
    [num1,txt1,raw1]=xlsread(xlsname,'materials');
    [num2,txt2,raw2]=xlsread(xlsname,'geometry');
    if any(strcmpi('bend & sweep',sheets))
        [~,txt3,raw3]=xlsread(xlsname,'bend & sweep');
        num3=cell2mat(raw3(4:end,1:6)); % ensure num3 has all six columns
    else
        num3=zeros(1,6);
        % specify default method and type
        txt3={'Pre-Bend','method ->','shear','Sweep','method ->','normal';
              '','type ->','poly','','type ->','poly'};
    end
catch ME
    fprintf('Could not read "%s". Please check format.\n',xlsname);
    rethrow(ME);
end
row1=4; % first row of useful data in the Excel file, materials sheet
row2=4; % first row of useful data in the Excel file, geometry sheet
row3=4; % first row of useful data in the Excel file, bend & sweep sheet
projname=txt2{1,1};

output_file='MatDBsi.txt';
if ~exist('output_file','var')  %print to command window if no output_file given
    fid1=1;
else
    fid1=fopen(output_file,'wt');   %try to open output_file for Writing in Text mode
    if (fid1 == -1)
        error('Could not open file "%s"\n',output_file);
        return
    end
end

output_file=sprintf('%s.nmd',projname);
if ~exist('output_file','var')  %print to command window if no output_file given
    fid2=1;
else
    fid2=fopen(output_file,'wt');   %try to open output_file for Writing in Text mode
    if (fid2 == -1)
        error('Could not open file "%s"\n',output_file);
        return
    end
end

% Indecies to columns
% materials sheet
mattypeid=3;
matnameid=2;
matrefid=17;
matthickid=4;
matexid=5;
matdenid=14;
matnuxyid=11;
Nmat=num1(end,1); % total number of materials
%geometry sheet
spanid=2; % index to column with span distances
afid=3; % index to column with airfoil shape names
typeid=4; % index to column with trailing edge type in Excel sheet
twid=5; % index to column with twist data in Excel sheet
chid=6; % index to column with chord data in Excel sheet
xoffid=7; % index to column with x-offset data in Excel sheet
acid=8;  % index to column with aerodynamic center data in Excel sheet
nlayid=10; % index to column with information about the first stack
nstacks=num2(1,8); % number of stacks defined in Excel sheet
nsegments=num2(2,8); % number of chordwise segments in Excel sheet
nsw=num2(1,6);  % number of shear webs in Excel sheet
% bend & sweep sheet
curvespan=1;
curvedisp=2;
curveslope=3;
sweepspan=4;
sweepdisp=5;
sweepslope=6;
% other important parameters
npanels=40;  % for interpolation of existing shapes to determine new shapes
Nsta=num2(end,1); % read the total number of stations from Excel sheet
mat_id=num2(2,(nlayid:nlayid+nstacks-1)); % array of material ids associated with each stack
comp_names=txt2(3,(nlayid+nstacks+1:nlayid+nstacks+1+nsegments-1)); % cell array of names for stack composites
mat_names=[];  % initialize the active list of materials
mat_names=['**UNSPECIFIED**'];

% Error checking
tmp=diff(num1(:,1));
if find(tmp~=1)
    error('Material ID''s in Materials Sheet are not in ascending order beginning with 1')
end
if num1(1,1)~=1
    error('Material ID''s in Materials Sheet are not in ascending order beginning with 1')
end


% Write isotropic and orthotropic material to material data file
% Perform operations for each material
try
    for i=1:Nmat
        fprintf(fid1,' <material>\n');
        type=lower(txt1{row1+i-1,mattypeid});
        switch type
            case 'isotropic'
                fprintf(fid1,'    <type>isotropic</type>\n');
                fprintf(fid1,'    <name>%s</name>\n',txt1{row1+i-1,matnameid});
                fprintf(fid1,'    <reference>%s</reference>\n',txt1{row1+i-1,matrefid});
                fprintf(fid1,'    <ex>%6.4e</ex>\n',num1(i,matexid)*1e6);
                fprintf(fid1,'    <dens>%4.0f</dens>\n',num1(i,matdenid));
                fprintf(fid1,'    <nuxy>%4.2f</nuxy>\n',num1(i,matnuxyid));
            case 'orthotropic'
                fprintf(fid1,'    <type>orthotropic</type>\n');
                fprintf(fid1,'    <name>%s</name>\n',txt1{row1+i-1,matnameid});
                fprintf(fid1,'    <reference>%s</reference>\n',txt1{row1+i-1,matrefid});
                fprintf(fid1,'    <ex>%6.4e</ex>\n',num1(i,matexid)*1e6);
                fprintf(fid1,'    <ey>%6.4e</ey>\n',num1(i,matexid+1)*1e6);
                fprintf(fid1,'    <ez>%6.4e</ez>\n',num1(i,matexid+1)*1e6);  % ez=ey
                fprintf(fid1,'    <gxy>%6.4e</gxy>\n',num1(i,matexid+3)*1e6);
                fprintf(fid1,'    <gyz>%6.4e</gyz>\n',num1(i,matexid+3)*1e6); % gxz=gyz=gxy
                fprintf(fid1,'    <gxz>%6.4e</gxz>\n',num1(i,matexid+3)*1e6); % gxz=gyz=gxy
                fprintf(fid1,'    <prxy>%4.2f</prxy>\n',num1(i,matnuxyid));
                fprintf(fid1,'    <pryz>%4.2f</pryz>\n',num1(i,matnuxyid)); % nuxz=nuyz=nuxy
                fprintf(fid1,'    <prxz>%4.2f</prxz>\n',num1(i,matnuxyid)); % nuxz=nuyz=nuxy
                 fprintf(fid1,'    <dens>%4.0f</dens>\n',num1(i,matdenid));
            otherwise
                error('Material %i: Material type not defined correctly in Excel file',i); break;
        end
        fprintf(fid1,' </material>\n');
    end
catch ME
    if fid1~=1, fclose(fid1); end
	if fid2~=1, fclose(fid2); end
    fprintf('Trouble writing information for Material %i',i);
    rethrow(ME);
    return;
end

% write header information to NMD file
fprintf(fid2,'<numad release="v2.0.1">\n');
fprintf(fid2,' <blade>\n');
fprintf(fid2,'   <rotation>cw</rotation>\n');
fprintf(fid2,'   <presweep>\n');
fprintf(fid2,'     <method>%s</method>\n',txt3{1,sweepslope});
fprintf(fid2,'     <table>\n');
sweeptable=num3(:,sweepspan:sweepslope);
sweeptable(isnan(sweeptable(:,1)),:)=[]; % ignore rows without span value
fprintf(fid2,' %g %g %g\n',transpose(sweeptable));
fprintf(fid2,'     </table>\n');
fprintf(fid2,'     <pptype>%s</pptype>\n',txt3{2,sweepslope});
fprintf(fid2,'   </presweep>\n');
fprintf(fid2,'   <precurve>\n');
fprintf(fid2,'     <method>%s</method>\n',txt3{1,curveslope});
fprintf(fid2,'     <table>\n');
curvetable=num3(:,curvespan:curveslope);
curvetable(isnan(curvetable(:,1)),:)=[]; % ignore rows without span value
fprintf(fid2,' %g %g %g\n',transpose(curvetable));
fprintf(fid2,'     </table>\n');
fprintf(fid2,'     <pptype>%s</pptype>\n',txt3{2,curveslope});
fprintf(fid2,'   </precurve>\n');

% write station information to NMD file
try
    for i=1:Nsta
        fprintf(fid2,'   <station>\n');
        afname=txt2{row2+i-1,afid};
        aftype=txt2{row2+i-1,typeid};
        if strcmpi(aftype,''),error('Trailing edge type not defined in Excel file at station #%i',i);end
        twist=num2(row2+i-1,twid);
        locationz=num2(row2+i-1,spanid);
        xoffset=num2(row2+i-1,xoffid);
        aerocent=num2(row2+i-1,acid);
        chord=num2(row2+i-1,chid);
        types={};
        if strcmpi(afname,'interp') % determine if this is a defined shape or a shape for interpolation
            % get interpolated airfoil information
            % find indecies to adjacent defined shapes
            ids=findAdjacentAF(i,txt2);
            % resample the given shapes and station data
            afs=[];
            chords=[];
            twists=[];
            spanlocs=[];
            acs=[];
            xoffsets=[];
            for j=1:length(ids)
                af_in=readAirfoil(['airfoils/' txt2{row2+ids(j)-1,afid} '.txt']);
                if strcmpi(aftype,'round')
                    af_out=resampleAirfoil(af_in.coords, npanels, 'cosine');
                else
                    af_out=resampleAirfoil(af_in.coords, npanels, 'half-cosine');
                end
                afs=[afs af_out];
                chords=[chords num2(row2-1+ids(j),chid)];
                twists=[twists num2(row2-1+ids(j),twid)];
                spanlocs=[spanlocs num2(row2-1+ids(j),spanid)];
                xoffsets=[xoffsets num2(row2-1+ids(j),xoffid)];
                acs=[acs num2(row2-1+ids(j),acid)];
            end
            [newaf,newchord] = interpAirfoil( afs, chords, spanlocs, locationz);
            if strcmpi(aftype,'flat')
                newaf=newaf(1:end,:);
            else
                newaf=newaf(2:end-1,:);
            end
            % For the following three numbers, only determine the interpolated
            %  value if it hasn't already been defined in the Excel sheet
            if isnan(twist)
                newtwist=interp1(spanlocs,twists,locationz,'linear');
            else
                newtwist=twist;
            end
            if isnan(xoffset)
                newxoffset=interp1(spanlocs,xoffsets,locationz,'linear');
            else
                newxoffset=xoffset;
            end
            if isnan(aerocent)
                newac=interp1(spanlocs,acs,locationz,'linear');
            else
                newac=aerocent;
            end
            newafname=['Interp_' sprintf('%06.0f',locationz*1000)];
            writeNuMADAirfoil(newaf,sprintf('Interpolated airfoil created on %s for blade project %s',date,projname),['airfoils/' newafname '.txt'])
            fprintf(fid2,'     <airfoilname>%s</airfoilname>\n',newafname);
            fprintf(fid2,'     <tetype>%s</tetype>\n',aftype);
            fprintf(fid2,'     <degreestwist>%5.2f</degreestwist>\n',newtwist);
            fprintf(fid2,'     <locationz>%5.2f</locationz>\n',locationz);
            fprintf(fid2,'     <xoffset>%6.4f</xoffset>\n',newxoffset);
            fprintf(fid2,'     <aerocenter>%6.4f</aerocenter>\n',newac);
            fprintf(fid2,'     <chord>%6.4f</chord>\n',newchord);
            fprintf(fid2,'     <coords>\n');
            for j=1:length(newaf)
                fprintf(fid2,' %9.7f %9.7f\n',[newaf(j,1) newaf(j,2)]);
            end
        else
            fprintf(fid2,'     <airfoilname>%s</airfoilname>\n',afname);
            fprintf(fid2,'     <tetype>%s</tetype>\n',aftype);
            fprintf(fid2,'     <degreestwist>%5.2f</degreestwist>\n',twist);
            fprintf(fid2,'     <locationz>%5.2f</locationz>\n',locationz);
            fprintf(fid2,'     <xoffset>%6.4f</xoffset>\n',xoffset);
            fprintf(fid2,'     <aerocenter>%6.4f</aerocenter>\n',aerocent);
            fprintf(fid2,'     <chord>%6.4f</chord>\n',chord);
            fprintf(fid2,'     <coords>\n');
            af=readAirfoil(['airfoils/' afname '.txt']);
            %af=resampleAirfoil(af.coords, npanels, 'half-cosine');
            for j=1:length(af.coords)
                fprintf(fid2,' %9.7f %9.7f\n',[af.coords(j,1) af.coords(j,2)]);
            end
        end
        fprintf(fid2,'     </coords>\n');
        fprintf(fid2,'     <delineationpoint>\n');
        
        % get number of layers in each stack at this station
        nlay=num2(row2+i-1,(nlayid:nlayid+nstacks-1));
        % get segment boundaries at this station
        seg_start=num2(row2+i-1,(nlayid+nstacks:nlayid+nstacks+nsegments));
        % get the list of dp types for each dp at this station
        dptype=txt2(row2+i-1,(nlayid+nstacks+nsegments+1:nlayid+nstacks+nsegments+1+nsegments));
        % get the list of included stacks for each segment at this station
        incl_stacks=txt2(row2+i-1,(nlayid+nstacks+nsegments+1+nsegments+1:nlayid+nstacks+nsegments+1+nsegments+1+nsegments-1));
        
        for j=1:nsegments  % write composite for each segment at this station
            a=sscanf(incl_stacks{j},'%u,');
            if ~isempty(a)
                fprintf(fid1,' <material>\n');
                fprintf(fid1,'    <type>composite</type>\n');
                compname=sprintf('%06.0f_%s',locationz*1000,comp_names{j});
                mat_names=[mat_names,{compname}];
                fprintf(fid1,'    <name>%s</name>\n',compname);
                fprintf(fid1,'    <reference>Reference text</reference>\n');
                fprintf(fid1,'    <thicknessType>Constant</thicknessType>\n');
                nuniquelayers=length(a);
                fprintf(fid1,'    <uniqueLayers>%u</uniqueLayers>\n',nuniquelayers);
                fprintf(fid1,'    <symmetryType>none</symmetryType>\n');
                for k=1:nuniquelayers
                    fprintf(fid1,'    <layer>\n');
                    matid=mat_id(a(k));
                    fprintf(fid1,'       <layerName>%s</layerName>\n',txt1{row1+matid-1,matnameid});
                    numlayer=nlay(a(k));
                    if isnan(numlayer) % if the number of layers has not been
                        % specified, interpolate based on values that have been
                        % given inboard and outboard
                        col=nlayid+a(k)-1;
                        adjacentvals=findAdjacentNum(col,i,num2);
                        numlayer=interp1(adjacentvals(:,1),adjacentvals(:,2),num2(row2+i-1,2));
                        numlayer=round(numlayer);   % No fractional layers allowed
                    end
                    layerthick=num1(matid,matthickid);
                    thick=layerthick*0.001;  % remember to convert to meters from mm
                    fprintf(fid1,'       <thicknessA>%9.7f</thicknessA>\n',thick);
                    fprintf(fid1,'       <thicknessB>%9.7f</thicknessB>\n',thick);
                    fprintf(fid1,'       <quantity>%d</quantity>\n',numlayer);
                    fprintf(fid1,'       <theta>0</theta>\n');
                    fprintf(fid1,'    </layer>\n');
                end
                fprintf(fid1,' </material>\n');
            end
        end
        
        % Here's where I can add capability to specify DP types for other than "SINGLE"
        for j=1:nsegments+1
            if seg_start(j)~=9999
                if isnan(seg_start(j)) % if the location of segment boundary has not been
                    % specified, interpolate based on values that have been
                    % given inboard and outboard
                    col=nlayid+nstacks+j-1;
                    adjacentvals=findAdjacentNum(col,i,num2);
                    newx=num2(row2+i-1,2);
                    seg_start(j)=interp1(adjacentvals(:,1),adjacentvals(:,2),newx);
                end
                switch lower(dptype{j})
                    case 'double'
                        fprintf(fid2,' %6.4f double\n',seg_start(j));
                    case 'flare'
                        fprintf(fid2,' %6.4f single **UNSPECIFIED**\n',seg_start(j));
                        warning('DP Type "flare" is not currently supported by xls2nmd.m');
                    case 'hourglass'
                        fprintf(fid2,' %6.4f single **UNSPECIFIED**\n',seg_start(j));
                        warning('DP Type "hourglass" is not currently supported by xls2nmd.m');
                    otherwise
                        fprintf(fid2,' %6.4f single **UNSPECIFIED**\n',seg_start(j));
                end
            end
        end
        
        fprintf(fid2,'     </delineationpoint>\n');
        fprintf(fid2,'     <surfacematerial>\n');
        if i~=Nsta
            for j=1:nsegments
                a=sscanf(incl_stacks{j},'%u,');
                if ~isempty(a)
                    compname=sprintf('%06.0f_%s',locationz*1000,comp_names{j});
                    fprintf(fid2,' %s\n',compname);
                end
            end
        end
        fprintf(fid2,'     </surfacematerial>\n');
        fprintf(fid2,'   </station>\n');
    end
catch ME
    if fid1~=1, fclose(fid1); end
    if fid2~=1, fclose(fid2); end
    fprintf('Trouble writing information for Station %i',i);
    rethrow(ME);
end


% Work on shear webs
try
    for i=1:Nsta
        for ii=1:nsw
            % get number of layers in each stack
            nlay=num2(row2+i-1,(nlayid:nlayid+nstacks-1));
            % get the list of included stacks
            incl_stacks=txt2(row2+i-1,nlayid+nstacks+nsegments+nsegments+nsegments+2+2*(ii-1));
            % get dp connections for the shear web
            dp_conn=txt2(row2+i-1,nlayid+nstacks+nsegments+nsegments+nsegments+3+2*(ii-1));
            locationz=num2(3+i,spanid);
            
            if ~strcmpi(incl_stacks{1},'')
                a=sscanf(incl_stacks{1},'%u,');
                b=sscanf(dp_conn{1},'%u,');
                if b(1)<=b(2),error('DP Connections for SW#%i are out of order',ii);end
                % write composite for this section of shearweb
                fprintf(fid1,' <material>\n');
                fprintf(fid1,'    <type>composite</type>\n');
                compname=sprintf('%06.0f_SW%i',locationz*1000,ii);
                mat_names=[mat_names,{compname}];
                fprintf(fid1,'    <name>%s</name>\n',compname);
                fprintf(fid1,'    <reference>Reference text</reference>\n');
                fprintf(fid1,'    <thicknessType>Constant</thicknessType>\n');
                nuniquelayers=length(a);
                fprintf(fid1,'    <uniqueLayers>%u</uniqueLayers>\n',nuniquelayers);
                fprintf(fid1,'    <symmetryType>none</symmetryType>\n');
                for k=1:nuniquelayers
                    fprintf(fid1,'    <layer>\n');
                    matid=mat_id(a(k));
                    fprintf(fid1,'       <layerName>%s</layerName>\n',txt1{row1+matid-1,matnameid});
                    numlayer=nlay(a(k));
                    if isnan(numlayer)
                        % interpolate based on known values and given span location
                        col=nlayid+a(k)-1;
                        adjacentvals=findAdjacentNum(col,i,num2);
                        numlayer=interp1(adjacentvals(:,1),adjacentvals(:,2),num2(row2+i-1,2));
                        numlayer=round(numlayer);
                    end
                    layerthick=num1(matid,matthickid);
                    thick=layerthick*0.001;
                    fprintf(fid1,'       <thicknessA>%9.7f</thicknessA>\n',thick);
                    fprintf(fid1,'       <thicknessB>%9.7f</thicknessB>\n',thick);
                    fprintf(fid1,'       <quantity>%d</quantity>\n',numlayer);
                    fprintf(fid1,'       <theta>0</theta>\n');
                    fprintf(fid1,'    </layer>\n');
                end
                fprintf(fid1,' </material>\n');
                
                fprintf(fid2,'   <shearweb>\n');
                fprintf(fid2,'     <material>%s</material>\n',compname);
                fprintf(fid2,'     <beginstation>%i</beginstation>\n',i);
                fprintf(fid2,'     <endstation>%i</endstation>\n',i+1);
                fprintf(fid2,'     <corner>%i %i %i %i</corner>\n',b);
                fprintf(fid2,'   </shearweb>\n');
            end
        end
    end
catch ME
    if fid1~=1, fclose(fid1); end
    if fid2~=1, fclose(fid2); end
    fprintf('Trouble writing information for Shear Web %i at Station %i',ii,i)
    rethrown(ME);
    return;
end

fprintf(fid2,' </blade>\n');
fprintf(fid2,' <activematerials>\n');
fprintf(fid2,'   <list>\n');
for i=1:length(mat_names)
    fprintf(fid2,' %s\n',mat_names{i});
end
fprintf(fid2,'   </list>\n');
fprintf(fid2,'   <colors>\n');

activecolors={[255   0   0]/255;  % red
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
kclr = rem((1:numel(mat_names))-1,numel(activecolors))+1;
colors = activecolors(kclr);

for i=1:length(mat_names)
    fprintf(fid2,' %g %g %g\n',colors{i});
end
fprintf(fid2,'   </colors>\n');
fprintf(fid2,' </activematerials>\n');
fprintf(fid2,' <ansys>\n');
fprintf(fid2,'   <boundarycondition>cantilevered</boundarycondition>\n');
fprintf(fid2,'   <elementsystem>shell181</elementsystem>\n');
fprintf(fid2,'   <meshing>elementsize</meshing>\n');
fprintf(fid2,'   <smartsize>5</smartsize>\n');
meshsize=max(num2(:,chid))*0.025; % set mesh size to 5% of max chord
fprintf(fid2,'   <elementsize>%6.4f</elementsize>\n',meshsize);
fprintf(fid2,'   <shell7gen>1</shell7gen>\n');
fprintf(fid2,'   <dbgen>1</dbgen>\n');
fprintf(fid2,' </ansys>\n');
fprintf(fid2,' </numad>\n');

if fid1~=1, fclose(fid1); end
if fid2~=1, fclose(fid2); end

if fid2~=1
%     nargoutchk(0,1)
    if isequal(1,nargout)
        varargout{1} = sprintf('%s.nmd',projname);
    end
end

end

% ================================================================================

function ids=findAdjacentAF(i,txt2)
row2=4;
N=size(txt2,1)-row2+1;
id=[];
for j=1:N
    if ~strcmpi(txt2{row2+j-1,3},'interp');id=[id j];end
end
outpointers=find(id>i);
inpointers=find(id<i);
if (length(outpointers)==1 || length(inpointers)==1)
    ids=[id(inpointers(end)) id(outpointers(1))];
else
    ids=[id(inpointers(end-1)) id(inpointers(end)) id(outpointers(1)) id(outpointers(2))];
end
end

function adjacentvals=findAdjacentNum(col,i,num2)
row2=4;
N=size(num2,1)-row2+1;
id=[];
for j=1:N
    if ~isnan(num2(row2+j-1,col));id=[id j];end
end
outpointers=find(id>i);
inpointers=find(id<i);
ids=[id(inpointers(end)) id(outpointers(1))];
adjacentvals=[num2(row2+ids(1)-1,2) num2(row2+ids(1)-1,col); ...
    num2(row2+ids(2)-1,2) num2(row2+ids(2)-1,col)];
end
