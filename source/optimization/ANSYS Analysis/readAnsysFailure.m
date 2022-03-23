function [elemFailure] = readAnsysFailure(fileName)
    fid = fopen(fileName);
    terminate = 0;
    elemFailure = [];
    while(terminate == 0)
        fLine = fgetl(fid);
        if(fLine == -1)
            terminate = 1;
        else
            fArray = str2num(fLine);
            if(length(fArray) > 2)
                elemFailure = [elemFailure;fArray(2)];
                %elemFailure = elemFailure + 60000*fArray(2)^2;
            end
        end
    end
    fclose(fid);

%     fid = fopen(fileName);
%     terminate = 0;
%     elemFailure = zeros(numEls,1);
%     k = 1;
%     while(terminate == 0)
%         fLine = fgetl(fid);
%         if(fLine == -1)
%             terminate = 1;
%         else
%             fArray = str2num(fLine);
%             if(length(fArray) > 2)
%                 elemFailure(k) = fArray(2);
%             end
%         end
%     end
%     fclose(fid);
end

