function a = loadAscii(filename)
    fid = fopen(filename);
    if fid==-1
       error([filename ': File does not exist']);
    end
    numlines = 0;
    try
       while(  1 )
           tline = fgetl(fid);
           if ~ischar(tline), break, end
           numlines=numlines+1;
    end
       fseek(fid, 0, 'bof');
       a = fscanf(fid,'%g',[numlines inf]);
       [r c] = size(a);
       % The next two steps are required since MATLAB reads the data in column major format
       % and fscanf reads data in row major format
       a = reshape(a,c,r);
       a = a';
       fclose(fid);
    catch
       fclose(fid);
       error([filename ': Error while reading file']);
    end
