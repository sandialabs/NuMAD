function varargout = readBModesOut(input_file,nrows)
%READBMODESOUT  Read a BModes output file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   bmodes = readBModesOut(input_file_name,nrows)
%     Compatible with BModes v3.00 blade file format
% 
%     Returns a structure containing the BModes output data.
%     nrows is the number of rows in each table

    %% open file
    fid=fopen(input_file);
    if (fid == -1)
        error('Could not open input "%s"\n',input_file);
        return
    end

    %% get file description from header
    line=fgetl(fid); 
    bmodes.header=fgetl(fid);
    for i=1:6
        line=fgetl(fid);
    end
    
    k=1;
    temp=zeros(nrows,6);
    while (1)
        line=fgetl(fid);
        bmodes.freq(k) = sscanf(line(33:end),'%g',inf);
        for i=1:3
            line=fgetl(fid);
        end
        for i=1:nrows
            line=fgetl(fid);
            temp(i,:)=sscanf(line,'%g',[1,6]); 
        end
        bmodes.tab{k}=temp;
        line=fgetl(fid);
        if (numel(line)~=0 && line(1)=='=');
            break;
        end
        line=fgetl(fid);
        k=k+1;
    end

    %% close file
    status=fclose(fid);

    if (nargout==0)
        numtab = numel(bmodes.tab);
        for k=1:numtab
            subplot(ceil(numtab/5),5,k);
            plot(bmodes.tab{k}(:,1),bmodes.tab{k}(:,2),...
                 bmodes.tab{k}(:,1),bmodes.tab{k}(:,4),...
                 bmodes.tab{k}(:,1),bmodes.tab{k}(:,6));
            if k==numtab, legend('flap','lag','twist'); end
        end
    else
        varargout{1} = bmodes;
    end
%endfunction