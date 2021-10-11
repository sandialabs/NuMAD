function data = read_1_ANSYSoutput(filename)

fid = fopen(filename,'r');
if (fid == -1)
    error('Could not open file "%s"',filename);
end


format='%s %f %s'; %Initialzie variable

file = fileread(filename); %Read whole file to determine max number of lines

lines=strsplit(file,'\n'); %Split file into individual lines (after every new line)
imax=length(lines);        %The max number of lines in the file.


data=inf; %Initialize variable to store numeric data. 

ct=0; %Counter for the numer of times data is written to tempdata
for i=1:imax  %For every line in the file
    tline = fgetl(fid); %read the next line
    try
       a=textscan(tline,format); % A try command is needed since textscan
                                 % will only work if "format" is readable
                                 % in that line. 
       if ~isempty(a{2})
           data=a{2};   %Overwrite any previous value of val so that only the last line with
                       % w/ a numerical argument is reported.
       end
    catch 
        %do nothing if textscan does did not work
    end
end
fclose(fid);