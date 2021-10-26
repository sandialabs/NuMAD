function varargout = readAirfoilData(input_file,varargin)
%READAIRFOILDATA  Read an AeroDyn airfoil file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   airfoil = readAirfoilData('input_file_name')
%   Returns a structure containing the data in an airfoil data file.
%
%   readAirfoilData('input_file_name',[AoA_limits])
%   Calling without an output argument will plot the data.
%   Specifying the 2x1 or 1x2 angle-of-attack limits is optional.


% open file
fid=fopen(input_file);
if (fid == -1)
    error('Could not open input "%s"\n',input_file);
    return
end

airfoil.title{1}=fgetl(fid);
airfoil.title{2}=fgetl(fid);
NTables=fscanf(fid,'%d',[1 1]); line=fgetl(fid);
airfoil.NTables=NTables;
airfoil.TableIDs=fscanf(fid,'%g',[1 NTables]); line=fgetl(fid);
line=fgetl(fid);  % unused input line
line=fgetl(fid);  % unused input line
line=fgetl(fid);  % unused input line
line=fgetl(fid);  % unused input line
airfoil.ALPHAL=fscanf(fid,'%g',[1 NTables]); line=fgetl(fid);
airfoil.CNA=fscanf(fid,'%g',[1 NTables]); line=fgetl(fid);
airfoil.CNS=fscanf(fid,'%g',[1 NTables]); line=fgetl(fid);
airfoil.CNSL=fscanf(fid,'%g',[1 NTables]); line=fgetl(fid);
airfoil.AOD=fscanf(fid,'%g',[1 NTables]); line=fgetl(fid);
airfoil.CDO=fscanf(fid,'%g',[1 NTables]); line=fgetl(fid);

i=1;
while (1)
    line=fgetl(fid);
    if (line == -1)
        break
    elseif isempty(line)
        break
    end
    [data,count]=sscanf(line,'%g',inf);
    if (count==0), break, end;
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

for n=1:NTables
    if (NCx==3)
        airfoil.table{n}=[AoA' CL(n,:)' CD(n,:)' CM(n,:)'];
    else
        airfoil.table{n}=[AoA' CL(n,:)' CD(n,:)'];
    end
end

airfoil.AoA=AoA;
airfoil.CL=CL;
airfoil.CD=CD;
if (NCx==3)
    airfoil.CM=CM;
end

%% close file
status=fclose(fid);

if (nargout==0)
    range=1:length(airfoil.AoA);
    if (nargin==2)
        AoA_limits=varargin{1};
        n1=find(airfoil.AoA>=AoA_limits(1),1);
        n2=find(airfoil.AoA>=AoA_limits(2),1);
        if isempty(n1), n1=range(1), end
        if isempty(n2), n2=range(2), end
        range = n1:n2;
    end
    subplot(NCx+1,1,1);
    if (NTables>4)
        set(gcf,'DefaultAxesColorOrder',prism(NTables));
    end
    plot(airfoil.AoA(range),airfoil.CL(:,range)); legend(num2str(airfoil.TableIDs','%g'));
    ylabel('C_L'); xlabel('Angle of Attack');
    subplot(NCx+1,1,2);
    plot(airfoil.AoA(range),airfoil.CD(:,range));
    ylabel('C_D'); xlabel('Angle of Attack');
    subplot(NCx+1,1,3);
    plot(airfoil.AoA(range),airfoil.CL(:,range) ./ airfoil.CD(:,range));
    ylabel('C_L/C_D'); xlabel('Angle of Attack');
    if (NCx==3)
        subplot(NCx+1,1,4);
        plot(airfoil.AoA(range),airfoil.CM(:,range));
        ylabel('C_M'); xlabel('Angle of Attack');
    end
else
    varargout{1} = airfoil;
end

%endfunction