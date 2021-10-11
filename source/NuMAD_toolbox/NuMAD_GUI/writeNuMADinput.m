function writeNuMADinput(data,filename)
%WRITENUMADINPUT  Output a NuMAD project file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeNuMADinput(data,filename)
%     Save NuMAD project data out to project file.

%===== CREDITS & CHANGELOG ================================================
%Developed by Wind & Water Power Technologies, Sandia National Laboratories
%2011.01.19  JCB: first draft

% Open the file for Writing in Text mode
fid = fopen(filename,'wt');
if (fid == -1)
    error('Could not open file "%s"',filename);
end

fprintf(fid,'<numad release="%s">\n',data.ReleaseVersion);

fprintf(fid,'<blade>\n');
fprintf(fid,'  <rotation>%s</rotation>\n',data.blade.BladeRotation);
fprintf(fid,'  <presweep>\n');
fprintf(fid,'    <method>%s</method>\n',data.blade.PresweepRef.method);
fprintf(fid,'    <table>\n');
fprintf(fid,'%g %g %g\n',data.blade.PresweepRef.table');
fprintf(fid,'    </table>\n');
fprintf(fid,'    <pptype>%s</pptype>\n',data.blade.PresweepRef.pptype);
fprintf(fid,'  </presweep>\n');
fprintf(fid,'  <precurve>\n');
fprintf(fid,'    <method>%s</method>\n',data.blade.PrecurveRef.method);
fprintf(fid,'    <table>\n');
fprintf(fid,'%g %g %g\n',data.blade.PrecurveRef.table');
fprintf(fid,'    </table>\n');
fprintf(fid,'    <pptype>%s</pptype>\n',data.blade.PrecurveRef.pptype);
fprintf(fid,'  </precurve>\n');

station = data.station;
for k = 1:numel(station)
    fprintf(fid,'  <station>\n');
    fprintf(fid,'    <airfoilname>%s</airfoilname>\n',station(k).AirfoilName);
    fprintf(fid,'    <tetype>%s</tetype>\n',station(k).TEtype);
    fprintf(fid,'    <degreestwist>%g</degreestwist>\n',station(k).DegreesTwist);
    fprintf(fid,'    <locationz>%g</locationz>\n',station(k).LocationZ);
    fprintf(fid,'    <xoffset>%g</xoffset>\n',station(k).Xoffset);
    fprintf(fid,'    <aerocenter>%g</aerocenter>\n',station(k).AeroCenter);
    fprintf(fid,'    <chord>%g</chord>\n',station(k).Chord);
    %
    fprintf(fid,'    <coords>\n');
    for j = 1:size(station(k).coords,1)
        fprintf(fid,'%g %g\n',station(k).coords(j,:));
    end
    fprintf(fid,'    </coords>\n');
    %
    fprintf(fid,'    <delineationpoint>\n');
    for j = 1:numel(station(k).dp)
        fprintf(fid,'%g %s %s\n',station(k).dp(j),station(k).dptype{j},station(k).dpmaterial{j});
    end
    fprintf(fid,'    </delineationpoint>\n');
    %
    fprintf(fid,'    <surfacematerial>\n');
    if k < numel(station)
        for j = 1:numel(station(k).sm)
            fprintf(fid,'%s\n',station(k).sm{j});
        end
    end
    fprintf(fid,'    </surfacematerial>\n');
    %
    fprintf(fid,'  </station>\n');
end

shearweb = data.shearweb;
for k = 1:numel(shearweb)
    fprintf(fid,'  <shearweb>\n');
    fprintf(fid,'    <material>%s</material>\n',shearweb(k).Material);
    fprintf(fid,'    <beginstation>%d</beginstation>\n',shearweb(k).BeginStation);
    fprintf(fid,'    <endstation>%d</endstation>\n',shearweb(k).EndStation);
    fprintf(fid,'    <corner>%d %d %d %d</corner>\n',shearweb(k).Corner);
    fprintf(fid,'  </shearweb>\n');
end

fprintf(fid,'</blade>\n');

active = data.active;
fprintf(fid,'<activematerials>\n');
fprintf(fid,'  <list>\n');
for k = 1:numel(active.list)
    fprintf(fid,'%s\n',active.list{k});
end
fprintf(fid,'  </list>\n');
% fprintf(fid,'  <colors>\n');
% for k = 1:numel(active.color)
%     fprintf(fid,'%g %g %g\n',active.color{k});
% end
% fprintf(fid,'  </colors>\n');
fprintf(fid,'</activematerials>\n');

ansys = data.ansys;
fprintf(fid,'<ansys>\n');
fprintf(fid,'  <boundarycondition>%s</boundarycondition>\n',ansys.BoundaryCondition);
fprintf(fid,'  <elementsystem>%s</elementsystem>\n',ansys.ElementSystem);
fprintf(fid,'  <multiplelayerbehavior>%s</multiplelayerbehavior>\n',ansys.MultipleLayerBehavior);
fprintf(fid,'  <meshing>%s</meshing>\n',ansys.meshing);
fprintf(fid,'  <smartsize>%d</smartsize>\n',ansys.smartsize);
fprintf(fid,'  <elementsize>%g</elementsize>\n',ansys.elementsize);
fprintf(fid,'  <shell7gen>%d</shell7gen>\n',ansys.shell7gen);
fprintf(fid,'  <dbgen>%d</dbgen>\n',ansys.dbgen);
if isfield(ansys,'FailureCriteria')
    fprintf(fid,'  <failurecriteria>\n');
    for kfc = 1:length(ansys.FailureCriteria)
        if ansys.FailureCriteria{kfc,2}
            fprintf(fid,'    %s\n',ansys.FailureCriteria{kfc,1});
        end
    end
    fprintf(fid,'  </failurecriteria>\n');
end
fprintf(fid,'</ansys>\n');

if isfield(data,'plot3d')
    plot3d = data.plot3d;
    fprintf(fid,'<plot3d>\n');
    fprintf(fid,'  <n_panels>%d</n_panels>\n',plot3d.n_panels);
    fprintf(fid,'  <spacing>%s</spacing>\n',plot3d.spacing);
    fprintf(fid,'  <interpmethod>%s</interpmethod>\n',plot3d.interpmethod);
    fprintf(fid,'  <breakpoints>\n');
    fprintf(fid,'    %g\n',plot3d.breakpoints);
    fprintf(fid,'  </breakpoints>\n');
    fprintf(fid,'  <newspanloc>\n');
    fprintf(fid,'    %g\n',plot3d.newspanloc);
    fprintf(fid,'  </newspanloc>\n');
    fprintf(fid,'</plot3d>\n');
end

if isfield(data,'flutterInput')
    flutterInput = data.flutterInput;
    fprintf(fid,'<flutterinput>\n');
    fprintf(fid,'  <fstfile>%s</fstfile>\n',flutterInput.fstFile);
    fprintf(fid,'  <bladefile>%s</bladefile>\n',flutterInput.bladeFile);
    fprintf(fid,'  <aerofile>%s</aerofile>\n',flutterInput.aeroFile);
    fprintf(fid,'  <outfile>%s</outfile>\n',flutterInput.outFile);
    fprintf(fid,'  <lcs>\n');
    fprintf(fid,'    %g\n',flutterInput.LCS);
    fprintf(fid,'  </lcs>\n');
    fprintf(fid,'  <omegastart>%g</omegastart>\n',flutterInput.OmegaStart);
    fprintf(fid,'  <omegaincrement>%g</omegaincrement>\n',flutterInput.OmegaIncrement);
    fprintf(fid,'  <omegaend>%g</omegaend>\n',flutterInput.OmegaEnd);
    fprintf(fid,'</flutterinput>\n');
end



fprintf(fid,'</numad>\n');

% Close the file
fclose(fid);
