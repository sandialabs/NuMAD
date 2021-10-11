function out = readWTPerf(filename)
%READWTPERF  Read WT_Perf output files.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   
fr = fileReader(filename); % create a file reader object
hdrInfo = parseHeader(fr.getLine); % get header info
% abort if result type not identified
assert(~isempty(hdrInfo.type),'Does not appear to be a WT_Perf results file:\n"%s"',filename);
assert(~strcmpi(hdrInfo.type,'unknown'),...
    'WT_Perf results files must begin with "Results" or "Blade-element":\n%s',filename);

fr.discardLine; % "Generated on ..."
fr.discardLine; % "Input file title:"
out.title = strtrim(fr.getLine);
fr.discardLine; % fifth line should be empty
line6 = fr.getLine;
fr.pushLine(line6); % push back so that all read functions start at line 6

WTPver = hdrInfo.verMajor;
switch hdrInfo.type
    case 'results'
        if isempty(strtrim(line6))
            out.type = 'parametric';
            out.parametric = readParametric(fr,WTPver);
        else
            out.type = 'cases';
            out.cases = readCases(fr,WTPver);
        end
    case 'bed'
        out.type = 'bed';
        out.bed = readBED(fr,WTPver);
end

end

function hdrInfo = parseHeader(hdr)

    % === Regular expression parts ===
    % From beginning of string: ^
    % Match+name everything before first parenthesis: (?<type>[^\(]*)
    % Match opening parenthesis: \(
    % Match+name everything that's not a comma: (?<ver>[^,]*)
    % Match everything that's not a double quote: [^"]*
    % Match first double quote: "
    % Match+name everything that's not a double quote: (?<name>[^"]*)
    % ================================
    n = regexp(hdr,'^(?<type>[^\(]*)\((?<ver>[^,]*)[^"]*"(?<name>[^"]*)','names');
    
    hdrInfo = struct('type','','verFull','','verMajor','','inputFile','');
    if isempty(n) || isempty(n.type)
        hdrInfo.type = '';
        return
    end
    

    if strfind(n.type,'WT_Perf')
        if strncmpi(n.type,'results',7)
            hdrInfo.type = 'results';
        elseif strncmpi(n.type,'blade-element',7)
            hdrInfo.type = 'bed';
        else
            hdrInfo.type = 'unknown';
        end
    else
        hdrInfo.type = '';
    end
    
    hdrInfo.verFull = n.ver;
    hdrInfo.verMajor = regexp(n.ver,'v[\d+\.]*','match','once');
    hdrInfo.inputFile = n.name;

end

function cases = readCases(fr,WTPver)
    varNames = fr.getLine;
    varUnits = fr.getLine;
    cases.VariableNames = regexp(varNames,'[^\t]+','match');
    cases.VariableUnits = regexp(varUnits,'[^\t]+','match');
    
    switch WTPver
        case 'v3.05.00'
            % discard the 'CavConverge' column
            formatSpec = '%f%f%f%f%f%f%f%f%f%f%*s%*s%*[^\n\r]';
            cases.VariableNames = cases.VariableNames(1:end-1);
            cases.VariableUnits = cases.VariableUnits(1:end-2);
        otherwise
            error('Version %s not yet supported.',WTPver);
    end

    % Initialize variables for textscan.
    delimiter = '\t';
    startRow = 1;
    endRow = inf;
    
    % Read columns of data according to format string.
    dataArray = textscan(fr.fileID, formatSpec, endRow(1)-startRow(1)+1, ...
        'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', true);
    
    % Convert from cell to array of doubles
    cases.table = [dataArray{1:end}];
end

function tables = readParametric(fr,WTPver)
    fr.discardLine; % line 6
    fr.discardLine; % line 7
    kTable = 1;
    notDone = true;
    while notDone
        lineSeparator = fr.getLine; % (======) or (------)
        if isequal(lineSeparator,-1)
            % end of file
            break;
        end
        if length(lineSeparator)<3
            continue;
        end
        switch lineSeparator(1:3)
            case '==='
                fr.discardLine; % (blank)
                fr.discardLine; % (Conditions leading to the maximum Cp:)
                fr.discardLine; % (TSR = xxxx)
                fr.discardLine; % (Pitch = xxxx)
                fr.discardLine; % (MaxCp = xxxx)
                fr.discardLine; % (blank)
            case '---'
                tables(kTable) = readParTable(fr);
                kTable = kTable + 1;
        end
    end
    
end

function parTable = readParTable(fr)
    sectionHeader = fr.getLine;
    n = regexp(sectionHeader,'(?<OutVar>.*) for (?<TabPar>[^=]*) = (?<TabVal>.*)','names');
    parTable.OutputVariable = strtrim(regexp(n.OutVar,'[^\(]*','match','once'));
    parTable.TableParameter = strtrim(n.TabPar);
    TabVal = sscanf(n.TabVal,'%f',[1 1]);
    fr.discardLine; % (blank)
    
    rowcol = fr.getLine;
    n = regexp(rowcol,'(?<RowPar>\S+)\s+(?<ColPar>\S+)\s+\((?<ColUnits>[^\)]*)','names');
    parTable.RowParameter = n.RowPar;
    parTable.ColumnParameter = n.ColPar;
    parTable.TableValues = TabVal;
    
    firstline = fr.getLine;
    n = regexp(firstline,'\((?<RowUnits>[^\)]*)\)(?<ColVal>.*)','names');
    parTable.ColumnValues = sscanf(n.ColVal,'%f',[1 inf]);
    nCols = numel(parTable.ColumnValues);
    
    % Initialize variables for textscan.
    formatSpec = [repmat('%f',1,nCols+1), '%*[^\n\r]'];
    delimiter = '\t';
    startRow = 1;
    endRow = inf;
    
    % Read columns of data according to format string.
    dataArray = textscan(fr.fileID, formatSpec, endRow(1)-startRow(1)+1, ...
        'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', true);
    
    parTable.RowValues = [dataArray{1}];
    parTable.Table = [dataArray{2:end}];
    
end

function tables = readBED(fr,WTPver)
    fr.discardLine; % line 6
    kTable = 1;
    notDone = true;
    while notDone
        header = fr.getLine; % (Blade-element data for:)
        if isequal(header,-1)
            % end of file
            break;
        end
        if ~isempty(strfind(header,'element data'))
            par1 = regexp(fr.getLine,'\s*(?<val>\S+)\s+(?<name>.*)','names');
            par2 = regexp(fr.getLine,'\s*(?<val>\S+)\s+(?<name>.*)','names');
            par3 = regexp(fr.getLine,'\s*(?<val>\S+)\s+(?<name>.*)','names');
            tables(kTable).ParameterNames{1} = par1.name;
            tables(kTable).ParameterNames{2} = par2.name;
            tables(kTable).ParameterNames{3} = par3.name;
            tables(kTable).ParameterValues(1) = sscanf(par1.val,'%f',[1 1]);
            tables(kTable).ParameterValues(2) = sscanf(par2.val,'%f',[1 1]);
            tables(kTable).ParameterValues(3) = sscanf(par3.val,'%f',[1 1]);
            fr.discardLine;
            varNames = fr.getLine;
            varUnits = fr.getLine;
            tables(kTable).VariableNames = regexp(varNames,'[^\t]+','match');
            tables(kTable).VariableUnits = regexp(varUnits,'[^\t]+','match');
            
            switch WTPver
                case 'v3.05.00'
                    % strings on columns 17 and 24
                    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%s%*[^\n\r]';
                otherwise
                    error('Version %s not yet supported.',WTPver);
            end
            
            % Initialize variables for textscan.
            delimiter = '\t';
            startRow = 1;
            endRow = inf;
            
            % Read columns of data according to format string.
            dataArray = textscan(fr.fileID, formatSpec, endRow(1)-startRow(1)+1, ...
                'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', true);
            
            switch WTPver
                case 'v3.05.00'
                    % convert from T/F to 1/0
                    dataArray{17} = double(cellfun(@(x) upper(x)=='T',dataArray{17}));
                    dataArray{24} = double(cellfun(@(x) upper(x)=='T',dataArray{24}));
            end
            
            % Convert from cell to array of doubles
            tables(kTable).Table = [dataArray{1:end}];
            
            % increment index
            kTable = kTable + 1;
        end
    end
end