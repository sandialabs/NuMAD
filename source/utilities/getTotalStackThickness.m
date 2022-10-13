function thickness=getTotalStackThickness(stationStack)
    nLayers=length(stationStack.plygroups);
    thickness=0;
    for iLayer=1:nLayers
        thickness=thickness+stationStack.plygroups(iLayer).nPlies*stationStack.plygroups(iLayer).thickness/1000;
    end
end