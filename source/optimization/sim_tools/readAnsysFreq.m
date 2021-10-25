function [Freq] = readAnsysFreq(fileName)
    fid = fopen(fileName);
    fLine = fgetl(fid);
    while(~contains(fLine,'1'))
        fLine = fgetl(fid);
    end
    terminate = 0;
    modeNum = 0;
    Freq = [];
    while(terminate == 0)
        if(fLine == -1)
            terminate = 1;
        else
            lnLst = str2num(fLine);
            if(length(lnLst) >= 2)
			    if(lnLst(2) > 0.0)
                    modeNum = modeNum + 1;
                    Freq = [Freq,lnLst(2)];
				end
            end
            fLine = fgetl(fid);
        end
    end
    fclose(fid);
end

