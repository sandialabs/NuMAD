function [stations,shearwebs,sortorder,deletelist] = sort_stations(stations,shearwebs)
    nStations = numel(stations);
    LocationZ = zeros(1,nStations);
    for k = 1:nStations
        LocationZ(k) = stations(k).LocationZ;
    end
    [~, sortorder] = sort(LocationZ);
    stations = stations(sortorder);
    
    deletelist = [];
    for ksw = 1:numel(shearwebs)
        sw = shearwebs(ksw);
        sw.BeginStation = find(sortorder==sw.BeginStation);
        sw.EndStation = find(sortorder==sw.EndStation);
        if ((sw.EndStation-sw.BeginStation) ~= 1)
            deletelist(end+1) = ksw;
        else
            shearwebs(ksw) = sw;
        end
    end
    if ~isempty(deletelist)
        shearwebs(deletelist) = [];
    end
end