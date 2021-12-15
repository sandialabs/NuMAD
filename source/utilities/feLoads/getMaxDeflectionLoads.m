function loadsTable = getMaxDeflectionLoads(blade,IEC,IECOutput) 
    bladeLength = blade.ispan(end);
    halfdz=bladeLength/40;
    r=(halfdz:2*halfdz:bladeLength)';
    sAr = [IECOutput.MaxOoPDefl];
    maxCase = 1;
    maxDef = 0.0;
    for i1 = 1:length(sAr)
        if(abs(sAr(i1).data) > maxDef)
            maxDef = abs(sAr(i1).data);
            maxCase = i1;
        end
    end
    loadsPointer.time = sAr(maxCase).time;
    loadsPointer.file = sAr(maxCase).File;
    loadsTable = getForceDistributionAtTime(loadsPointer,IEC,r);
end

