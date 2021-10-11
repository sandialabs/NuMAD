function ConvertAD13to12(input_file,output_file,UseCM)
%CONVERTAD13TO12  Convert an AeroDyn v13.0 airfoil file to v12.5
%                          Under Construction
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   ConvertAD13to12('input_file_name','output_file_name',UseCM)
%   Writes a v12.5 airfoil file using data and comments from a v13.0 file
%
%   UseCM==true  will include CM data (if available)
%   UseCM==false will not include CM data
%
%   Input file format defined in AeroDyn User's Guide, Version 12.50

%% read airfoil file
fid=fopen(input_file);
if (fid == -1)
    error('Could not open input "%s"',input_file);
    return
end

title{1}=fgetl(fid);
title{2}=fgetl(fid);
title{3}=fgetl(fid);
NTables=fscanf(fid,'%d',[1 1]); comment{1}=strtrim(fgetl(fid));  % Number of airfoil tables

for n=1:NTables
    MulTabMet(n)=fscanf(fid,'%g',[1 NTables]); comment{2}=strtrim(fgetl(fid));  % Table ID parameters
    ControlSetting = fgetl(fid);
    Stall(n)=fscanf(fid,'%g',[1 NTables]); comment{3}=strtrim(fgetl(fid));  % stall angle
    ALPHAL(n)=fscanf(fid,'%g',[1 NTables]); comment{4}=strtrim(fgetl(fid));  % zero-lift AoA
    CNA(n)=fscanf(fid,'%g',[1 NTables]); comment{5}=strtrim(fgetl(fid));  % Cn slope near zero lift
    CNS(n)=fscanf(fid,'%g',[1 NTables]); comment{6}=strtrim(fgetl(fid));  % Cn at positive static stall
    CNSL(n)=fscanf(fid,'%g',[1 NTables]); comment{7}=strtrim(fgetl(fid));  % Cn at negative static stall
    AOD(n)=fscanf(fid,'%g',[1 NTables]); comment{8}=strtrim(fgetl(fid));  % AoA for minimum drag
    CDO(n)=fscanf(fid,'%g',[1 NTables]); comment{9}=strtrim(fgetl(fid));  % minimum Cd
    
    i=1;
    while (1)
        line=fgetl(fid);
        if strfind(line,'EOT')
            break
        elseif (line == -1)
            status=fclose(fid);
            error('EOT not found.  End of file found instead');
        elseif isempty(line)
            status=fclose(fid);
            error('EOT not found.  Empty line found instead');
        end
        [data,count]=sscanf(line,'%g',inf);
        AoA(n,i)=data(1);
        CL(n,i)=data(2);
        CD(n,i)=data(3);
        if (count==4)
            CM(n,i)=data(4);
        end
        i=i+1;
    end
end %(for n=1:NTables)

status=fclose(fid);

for n=1:NTables
    if ~isequal(AoA(1,:),AoA(n,:))
        error('Each table must have the same angle of attack values');
    end
end

%% generate new airfoil file
% the 't' in 'wt' below has something to do with the way newlines
% are handled - I had to add this so that Notepad can read the file
% properly
fid=fopen(output_file,'wt');
if (fid == -1)
    error('Could not open output "%s"',output_file);
end
fprintf(fid,'%s\n',title{1});
fprintf(fid,'%s\n',title{2});
fprintf(fid,'%4d  %s\n',NTables,comment{1});
str = repmat([' %8g'],1,NTables);
fprintf(fid,str,MulTabMet);
fprintf(fid,'  %s\n',comment{2});
fprintf(fid,str,Stall);
fprintf(fid,'  %s\n',comment{3});
fprintf(fid,str,repmat([0],1,NTables));
fprintf(fid,'  No longer used, enter zero\n');
fprintf(fid,str,repmat([0],1,NTables));
fprintf(fid,'  No longer used, enter zero\n');
fprintf(fid,str,repmat([0],1,NTables));
fprintf(fid,'  No longer used, enter zero\n');
fprintf(fid,str,ALPHAL);
fprintf(fid,'  %s\n',comment{4});
fprintf(fid,str,CNA);
fprintf(fid,'  %s\n',comment{5});
fprintf(fid,str,CNS);
fprintf(fid,'  %s\n',comment{6});
fprintf(fid,str,CNSL);
fprintf(fid,'  %s\n',comment{7});
fprintf(fid,str,AOD);
fprintf(fid,'  %s\n',comment{8});
fprintf(fid,str,CDO);
fprintf(fid,'  %s\n',comment{9});

for i=1:length(AoA)
    fprintf(fid,'%7.2f',AoA(1,i));
    for n=1:NTables
        fprintf(fid,'   %6.3f %7.4f',CL(n,i),CD(n,i));
        if count==4 & UseCM
            fprintf(fid,' %7.4f',CM(n,i));
        end
    end
    if (i ~= length(AoA))
        fprintf(fid,'\n');
    end
end

status=fclose(fid);

%endfunction