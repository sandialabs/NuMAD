function loadNuMAD(arg)
% LOADNUMAD  Convert SData3.nmd and SData1.nmd into Matlab data structure.
%
%  Important Notes/Limitations:
%   -?
%

% ========================================================
%   Written by Brian Resor, Sandia National Laboratories
%   Last update: 01/21/2011

% get blade rotation direction from SData1.nmd
fid=fopen('SData1.nmd');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if strfind(tline,'set BladeRotation')
        tmp=regexp(tline,'[c+w+]*','match');
        tmp=tmp{1};
        switch tmp
            case 'cw'
                rotationdirection=-1;
            case 'ccw'
                rotationdirection=1;
            otherwise 
                    error('Rotation direction was not read properly')
        end
    end
end
fclose(fid);

fid=fopen('SData3.nmd');

while 1
    tline = fgetl(fid);
     if ~ischar(tline), break, end
    %% Read Shear Web information
    if strfind(tline,'set SW(')
        pat='\((?<number>\w+),(?<name>\w+)\)\s+(?<setting>\w+)|\((?<number>\w+),(?<name>\w+),(?<number2>\w+)\)\s+(?<setting>\w+)';
        n=regexp(tline,pat,'names');
        if ~isempty(n.number2)
            eval(['SW(' n.number ').' n.name '(' n.number2 ')=' n.setting ';']);
        else
            try
                eval(['SW(' n.number ').' n.name '=' n.setting ';']);
            catch
                eval(['SW(' n.number ').' n.name '=''' n.setting ''';']);
            end
        end
        
        %% Read SMN information
    elseif strfind(tline,'set SMN(')
        pat='\((?<number>\w+),(?<number2>\w+)\)\s+(?<setting>\w+)';
        n=regexp(tline,pat,'names');
        eval(['SMN{' num2str(n.number) '}.Matl{' n.number2 '}=''' n.setting ''';']);
        
        %% Read BSP information
    elseif strfind(tline,'set BSP(')
        pat='\((?<number>\w+),(?<number2>\w+),(?<name>\w+)\)\s+(?<setting>[-.0-9e]+)|\((?<number>\w+),(?<number2>\w+),(?<name>\w+)\)\s+(?<setting>\w+)|\((?<number>\w+),(?<name>\w+)\)\s+(?<setting>\w+)';
        n=regexp(tline,pat,'names');
        if ~isempty(n.number2)
            try
                eval(['StationData(' n.number ').' n.name '(' n.number2 ')=' n.setting ';']);
            catch
                eval(['StationData(' n.number ').' n.name '{' n.number2 '}=''' n.setting ''';']);
            end
        else
            eval(['StationData(' n.number ').' n.name '=' n.setting ';']);
        end
        
    elseif isempty(strfind(tline,'AirfoilName'))
        pat='(?<name>\w+)\((?<number>\w+)\)\s+(?<setting>[-.0-9e]+)|(?<name>\w+)\((?<number>\w+)\)\s+(?<setting>\w+)';
    n=regexp(tline,pat,'names');
        try
            eval(['StationData(' n.number ').' n.name '=' n.setting ';']);
        catch
            eval(['StationData(' n.number ').' n.name '=''' n.setting ''';']);
        end
    end
end

for i=1:length(StationData)
    StationData(i).rotdir=rotationdirection;
end
 
fclose(fid);

if exist('SW','var')
    save SData SW SMN StationData
else
    save SData SMN StationData
end

end