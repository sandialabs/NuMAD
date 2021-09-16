function write_plot3d(app,blade,filename)
%WRITE_PLOT3D  Exports blade surface geometry to Plot3D format. 
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   write_plot3d(app,blade,filename)
%     app = app data structure in NuMAD
%     blade = blade data structure in NuMAD
%     filename = plot3d file to be written
%   
%   Plot3D format options are stored in app.plot3d

%jcb: passing both "app" and "blade" is temporary until the data
%     re-organization is finalized

if isfield(app,'plot3d')
    n_panels = app.plot3d.n_panels;
    spacing = app.plot3d.spacing;
    interpmethod = app.plot3d.interpmethod;
    newspanloc = app.plot3d.newspanloc;
    breakpoints = app.plot3d.breakpoints;
else
    n_panels = 100;
    spacing = 'cosine';
    interpmethod = 'linear';
    newspanloc = []; 
    breakpoints = [];
end

twistFlag = 1;  % ccw rotor rotation
if strcmp(app.BladeRotation,'cw')
    twistFlag = -1;  % cw rotor rotation
end

% resample the airfoils with the specified number of panels
ResampledStation = resampleAirfoilDB(app.station, n_panels, spacing);
chordspacing = ResampledStation(1).c;
ResampledStation = rmfield(ResampledStation,...
    {'AirfoilName','TEtype','AeroCenter','coords','dp','dptype','dpmaterial','sm','LE','s','c'});
nResampledPoints = size(ResampledStation(1).xy,1);
TotalStations = numel(ResampledStation);

% create intermediate stations by interpolation
[spanloc,chord,xoffset,degreestwist] = deal(zeros(1,TotalStations));
for kStation = 1:TotalStations
    spanloc(kStation) = ResampledStation(kStation).LocationZ;
    chord(kStation) = ResampledStation(kStation).Chord;
    xoffset(kStation) = ResampledStation(kStation).Xoffset;
    degreestwist(kStation) = ResampledStation(kStation).DegreesTwist;
end
newstation = ResampledStation(1);  % ensure similar structures
for kNewSpan = 1:numel(newspanloc)
    kb = find(spanloc>newspanloc(kNewSpan),1);
    if isequal('spline',interpmethod) &&  kb>2 && kb<TotalStations-1
        ka = kb-2;
        kb = kb+1;
        newstation.Xoffset = interp1(spanloc,xoffset,newspanloc(kNewSpan),'spline');
        newstation.DegreesTwist = interp1(spanloc,degreestwist,newspanloc(kNewSpan),'spline');
    else
        ka = kb-1;
        newstation.Xoffset = interp1(spanloc,xoffset,newspanloc(kNewSpan),'linear');
        newstation.DegreesTwist = interp1(spanloc,degreestwist,newspanloc(kNewSpan),'linear');
    end
    %jcb: need to leave off last point to work with interpAirfoil
    coords = zeros(nResampledPoints-1,2*(kb-ka+1));
    j=1;
    for k=ka:kb
        coords(:,j:j+1) = ResampledStation(k).xy(1:end-1,:);
        j=j+2;
    end
    [newaf,newchord] = interpAirfoil(coords,chord(ka:kb),spanloc(ka:kb),newspanloc(kNewSpan));
    newstation.xy = [newaf; newaf(end,:)]; %jcb: add last point back on
    newstation.Chord = newchord;
    newstation.LocationZ = newspanloc(kNewSpan);
    ResampledStation(end+1) = newstation; %#ok
end

% re-sort the new stations
TotalStations = numel(ResampledStation);
spanloc = zeros(1,TotalStations);
for kStation = 1:TotalStations
    spanloc(kStation) = ResampledStation(kStation).LocationZ;
end
[~,sortvector] = sort(spanloc);
ResampledStation = ResampledStation(sortvector);

% initialize matrices for transformed x,y,z
[xt,yt,zt] = deal(zeros(nResampledPoints-2,TotalStations));
for kStation = 1:TotalStations
    station = ResampledStation(kStation);
    % jcb: need to add the flexibility of choosing resampling in the next
    % step
    xy = station.xy; % use resampled points
    % jcb: when using the xy data from NuMAD, the first and last 
    %      points are duplicated; therefore, need to start at 2 
    %      and stop at nResampledPoints-1
    x = (xy(2:nResampledPoints-1,1) - station.Xoffset) * station.Chord * twistFlag;
    y = (xy(2:nResampledPoints-1,2)                  ) * station.Chord;
    twist = twistFlag * station.DegreesTwist * pi/180;
    coords = zeros(length(x),4);
    coords(:,1) = cos(twist) * x - sin(twist) * y;
    coords(:,2) = sin(twist) * x + cos(twist) * y;
    %coords(:,3) = zeros(size(x));
    coords(:,4) = ones(size(x));
    %         xt = cos(twist) * x - sin(twist) * y;
    %         yt = sin(twist) * x + cos(twist) * y;
    %         zt = ones(size(x)) * station.LocationZ;
    %         xn = station.xn;
    
    % use the generating line to translate and rotate the coordinates
    presweep_slope = ppval(blade.PresweepRef.dpp,station.LocationZ);
    precurve_slope = ppval(blade.PrecurveRef.dpp,station.LocationZ);
    [presweep_rot, precurve_rot] = deal(0);
    if isequal(blade.PresweepRef.method,'normal')
        presweep_rot = atan(presweep_slope*twistFlag);
    end
    if isequal(blade.PrecurveRef.method,'normal')
        precurve_rot = atan(-precurve_slope);
    end
    transX = twistFlag*ppval(blade.PresweepRef.pp,station.LocationZ);
    transY = ppval(blade.PrecurveRef.pp,station.LocationZ);
    R = makehgtform('yrotate',presweep_rot,'xrotate',precurve_rot);
    T = makehgtform('translate',transX,transY,station.LocationZ);
    coords = coords * R' * T';
    % jcb: important for each station's data to be in a column 
    %      rather than a row; see fprintf statement below
    xt(:,kStation) = coords(:,1);
    yt(:,kStation) = coords(:,2);
    zt(:,kStation) = coords(:,3);
end

if (0)
    %%% DEBUGGING
    assignin('base','xt',xt); %#ok<UNRCH>
    assignin('base','yt',yt);
    assignin('base','zt',zt);
    mesh(zt,xt,yt); axis equal;
end

indicesOfBreakpoints = zeros(1,numel(breakpoints));
for kBreakpoint = 1:numel(breakpoints)
    bp = breakpoints(kBreakpoint);
    [~,ind] = min(sqrt((chordspacing-bp).^2));
    indicesOfBreakpoints(kBreakpoint) = ind - 1;  % see lines 91-92 above
end
indicesOfBreakpoints = unique([1, indicesOfBreakpoints, size(xt,1)]);

fid = fopen(filename,'wt'); % open for Writing in Text mode
if (fid == -1)
    errordlg(sprintf('Could not open file "%s"',filename),'error');
    error('Could not open file "%s"',filename);
end

% output the data in Plot3d format
try
    nBlocks = numel(indicesOfBreakpoints)-1; % number of xyz blocks
    fprintf(fid,'%d\n',nBlocks);
    for kblock = 1:nBlocks
        a = indicesOfBreakpoints(kblock);
        b = indicesOfBreakpoints(kblock+1);
        fprintf(fid,'%d  %d  %d\n',1+b-a,size(xt,2),1);
    end
    
    columnsPerLine = 5;
    for kblock = 1:nBlocks
        a = indicesOfBreakpoints(kblock);
        b = indicesOfBreakpoints(kblock+1);
        fprintf_matrix(fid,xt(a:b,:),columnsPerLine);
        fprintf_matrix(fid,yt(a:b,:),columnsPerLine);
        fprintf_matrix(fid,zt(a:b,:),columnsPerLine);  
    end
    
catch ME
    % The try..catch..end ensures the file gets closed
    % in case of a programming error.
    fclose(fid);
    rethrow(ME);
end
fclose(fid);
msgbox(sprintf('Plot3D output written to %s',filename),'Notification');

end % function

function fprintf_matrix(fid,matrixData,columnsPerLine)
    kColumn = 1;
    for kData=1:numel(matrixData)
        fprintf(fid,'%g ',matrixData(kData));
        kColumn = kColumn + 1;
        if kColumn > columnsPerLine
            fprintf(fid,'\n');
            kColumn = 1;
        end
    end
end
