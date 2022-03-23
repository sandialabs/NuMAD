function data = readANSYSoutputs(filename,ncol)

    fid = fopen(filename,'r');
    if (fid == -1)
        error('Could not open file "%s"',filename);
    end

 
    format=''; %Initialzie variable
    for i=1:ncol 
        format = strcat(format,'%f'); %The format for the number of columns
    end

    file = fileread(filename); %Read whole file to determine max number of lines

    lines=strsplit(file,'\n'); %Split file into individual lines (after every new line)
    imax=length(lines);        %The max number of lines in the file.

 
    tempdata=zeros(imax,ncol); %Matrix to store numeric data. After all lines 
    %have been read, then the number of lines with numeric data will be known. 
    %The remaining zeros will be truncated.

    ct=0; %Counter for the numer of times data is written to tempdata
    for i=1:imax  %For every line in the file
        tline = fgetl(fid); %read the next line
        try
            a=textscan(tline,format); % A try command is needed since textscan 
            b=cell2mat(a);                   % will only work if jmax floats are
                                        % in that line. 
            if ~isempty(b)
                ct=ct+1;
                tempdata(ct,:)=b;
            end
        catch 
            %do nothing if textscan does did not work
        end
%     disp(tline)
%     if ischar(tline)
%         ct=ct+1
%     end
    %data=[data; textscan(tline,'%f %f %f %f')];
    end
fclose(fid);
data=tempdata(1:ct,:); %Transver the read numeric data to a new variable.
end
