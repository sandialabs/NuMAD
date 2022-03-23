function linearLoadFactors = ReadANSYS_LinearBucklingResults(blade, config, iLoad, fid, bucklingFilename)

    fid=fopen([bucklingFilename '.out']);
    for jj=1:5
        tline = fgetl(fid);
    end
    data=cell(1,5);
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), 
            break, 
        end
        data=[data; textscan(tline,'%f %f %f %f %f')];
    end
    fclose(fid);
    disp(' ')
    data=cell2mat(data);
    linearLoadFactors=data(1:config.ansys.analysisFlags.globalBuckling,2); %Extract the load factors (LF) from the linear buckling analysis
    delete([bucklingFilename '.out'])
    
end
    