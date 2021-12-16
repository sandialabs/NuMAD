function writeMatDB(matdb,filename)
%WRITEMATDB  Output a NuMAD material database file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeMatDB(matdb,filename)
%     Output NuMAD material database FILENAME (typically 'MatDBsi.txt')
%

%===== CREDITS & CHANGELOG ================================================
%Developed by Wind & Water Power Technologies, Sandia National Laboratories
%2011.01.19  JCB: first draft based on writeNuMADinput.m

% Open the file for Writing in Text mode
fid = fopen(filename,'wt');
if (fid == -1)
    error('Could not open file "%s"',filename);
end

try
for k = 1:numel(matdb)
    fprintf(fid,'<material>\n');
    fprintf(fid,'   <type>%s</type>\n',matdb(k).type);
    fprintf(fid,'   <name>%s</name>\n',matdb(k).name);
    fprintf(fid,'   <reference>%s</reference>\n',matdb(k).reference);
    switch matdb(k).type
        case 'isotropic'
            fprintf(fid,'   <ex>%g</ex>\n',matdb(k).ex);
            fprintf(fid,'   <dens>%g</dens>\n',matdb(k).dens);
            fprintf(fid,'   <nuxy>%g</nuxy>\n',matdb(k).nuxy);
        case 'orthotropic'
            fprintf(fid,'   <ex>%g</ex>\n',matdb(k).ex);
            fprintf(fid,'   <ey>%g</ey>\n',matdb(k).ey);
            fprintf(fid,'   <ez>%g</ez>\n',matdb(k).ez);
            fprintf(fid,'   <gxy>%g</gxy>\n',matdb(k).gxy);
            fprintf(fid,'   <gyz>%g</gyz>\n',matdb(k).gyz);
            fprintf(fid,'   <gxz>%g</gxz>\n',matdb(k).gxz);
            fprintf(fid,'   <prxy>%g</prxy>\n',matdb(k).prxy);
            fprintf(fid,'   <pryz>%g</pryz>\n',matdb(k).pryz);
            fprintf(fid,'   <prxz>%g</prxz>\n',matdb(k).prxz);
            fprintf(fid,'   <dens>%g</dens>\n',matdb(k).dens);
        case 'composite'
            fprintf(fid,'   <thicknessType>%s</thicknessType>\n',matdb(k).thicknessType);
            fprintf(fid,'   <uniqueLayers>%d</uniqueLayers>\n',matdb(k).uniqueLayers);
            fprintf(fid,'   <symmetryType>%s</symmetryType>\n',matdb(k).symmetryType);
            for j = 1:matdb(k).uniqueLayers
                fprintf(fid,'   <layer>\n');
                fprintf(fid,'      <layerName>%s</layerName>\n',matdb(k).layer(j).layerName);
                fprintf(fid,'      <thicknessA>%g</thicknessA>\n',matdb(k).layer(j).thicknessA);
                fprintf(fid,'      <thicknessB>%g</thicknessB>\n',matdb(k).layer(j).thicknessB);
                fprintf(fid,'      <quantity>%g</quantity>\n',matdb(k).layer(j).quantity);
                fprintf(fid,'      <theta>%g</theta>\n',matdb(k).layer(j).theta);
                fprintf(fid,'   </layer>\n');
            end
        otherwise
            warning('Unknown material type: %s',matdb(k).type);
    end
    
    fcfields = {'xten','xcmp','yten','ycmp','zten','zcmp','xy','yz','xz',...
        'xycp','yzcp','xzcp','xzit','xzic','yzit','yzic','g1g2',...
        'etal','etat','alp0'};
    
    switch matdb(k).type
        case {'isotropic','orthotropic'}
            for kf = 1:length(fcfields)
                fcname = fcfields{kf};
                fcvalue = matdb(k).(fcname);
                if ~isempty(fcvalue)
                	fprintf(fid,'   <%s>%g</%s>\n',fcname,fcvalue,fcname);  
                %else
                %    fprintf(fid,'   <%s></%s>\n',fcname,fcname);  
                end
            end
    end

    fprintf(fid,'</material>\n');
end
% Close the file
fclose(fid);

catch ME
% Close the file
fclose(fid);
rethrow(ME);
end
end

