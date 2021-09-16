function ConvertAirfoilVersion(input_file,output_file)
%CONVERTAIRFOILVERSION  Convert an AeroDyn v12.5 airfoil file to v13.0
%                          Under construction
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   ConvertAirfoilVersion('input_file_name','output_file_name')
%   Writes a v13.0 airfoil file using data and comments from a v12.5 file
%
%   Input file format defined in AeroDyn User's Guide, Version 12.50

%% read airfoil file
fid=fopen(input_file);
if (fid == -1)
    fprintf('Error: could not open input "%s"\n',input_file);
    return
end

title{1}=fgetl(fid);
title{2}=fgetl(fid);
NTables=fscanf(fid,'%d',[1 1]); comment{1}=strtrim(fgetl(fid));  % Number of airfoil tables
MulTabMet=fscanf(fid,'%g',[1 NTables]); comment{2}=strtrim(fgetl(fid));  % Table ID parameters
Stall=fscanf(fid,'%g',[1 NTables]); comment{3}=strtrim(fgetl(fid));  % stall angle

for i=1:3
    line=fgetl(fid);  % unused input lines
end

ALPHAL=fscanf(fid,'%g',[1 NTables]); comment{4}=strtrim(fgetl(fid));  % zero-lift AoA
CNA=fscanf(fid,'%g',[1 NTables]); comment{5}=strtrim(fgetl(fid));  % Cn slope near zero lift
CNS=fscanf(fid,'%g',[1 NTables]); comment{6}=strtrim(fgetl(fid));  % Cn at positive static stall
CNSL=fscanf(fid,'%g',[1 NTables]); comment{7}=strtrim(fgetl(fid));  % Cn at negative static stall
AOD=fscanf(fid,'%g',[1 NTables]); comment{8}=strtrim(fgetl(fid));  % AoA for minimum drag
CDO=fscanf(fid,'%g',[1 NTables]); comment{9}=strtrim(fgetl(fid));  % minimum Cd

i=1;
while (1)
    line=fgetl(fid);
    if (line == -1)
        break
    elseif isempty(line)
        break
    end
    [data,count]=sscanf(line,'%g',inf);
    AoA(i)=data(1);
    NCx=(count-1)/NTables;  % equals 2 if data is given for CL and CD
    % equals 3 if data is given for CL, CD, and CM
    for n=1:NTables
        CL(n,i)=data(2+(n-1)*NCx);
        CD(n,i)=data(3+(n-1)*NCx);
        if (NCx==3)
            CM(n,i)=data(4+(n-1)*NCx);
        end
    end
    i=i+1;
end

status=fclose(fid);

title0='AeroDyn airfoil file.  Compatible with AeroDyn v13.0.';

%% generate new airfoil file
fid=fopen(output_file,'w');
if (fid == -1)
    fprintf('Error: could not open output "%s"\n',output_file);
    return
end
fprintf(fid,'%s\n',title0);
fprintf(fid,'%s\n',title{1});
fprintf(fid,'%s\n',title{2});
fprintf(fid,'%d   %s\n',NTables,comment{1});
for n=1:NTables
    fprintf(fid,' %8g  %s\n',MulTabMet(n),comment{2});
    fprintf(fid,' %8g  %s\n',Stall(n),comment{3});
    fprintf(fid,' %8g  %s\n',ALPHAL(n),comment{4});
    fprintf(fid,' %8g  %s\n',CNA(n),comment{5});
    fprintf(fid,' %8g  %s\n',CNS(n),comment{6});
    fprintf(fid,' %8g  %s\n',CNSL(n),comment{7});
    fprintf(fid,' %8g  %s\n',AOD(n),comment{8});
    fprintf(fid,' %8g  %s\n',CDO(n),comment{9});
    for i=1:length(AoA)
        fprintf(fid,' %4g %8g %8g',AoA(i),CL(n,i),CD(n,i));
        if (NCx==3)
            fprintf(fid,' %8g\n',CM(n,i));
        else
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'EOT\n');
end

status=fclose(fid);

%endfunction