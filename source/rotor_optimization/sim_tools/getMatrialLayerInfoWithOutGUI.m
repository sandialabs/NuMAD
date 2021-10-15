function [isoorthoInModel,compsInModel,SkinAreas,app] = getMatrialLayerInfoWithOutGUI(blade)
   %Temparary workaround to extract data needed for: 
    
    app=getApp('numad.nmd','ansys',blade.mesh);


    %From write_shell7.m
    TotalStations = numel(app.station);
    TotalShearwebs = numel(app.shearweb);
    SkinAreas(1:TotalStations-1) = struct('startIB',[],'endIB',[],'startOB',[],'endOB',[],'Material',{''});
    for kStation = 1:numel(SkinAreas)
        stationIB = app.station(kStation);
        stationOB = app.station(kStation+1);

        for kdp = 1:numel(stationIB.dp)-1
            switch stationIB.dptype{kdp}
                case {'single','double'}
                    % single and double are equivalent on the area inboard edge
                    % start and end are current and next DP
                    SkinAreas(kStation).startIB(end+1) = kdp;
                    SkinAreas(kStation).endIB(end+1)   = kdp+1;
                    SkinAreas(kStation).Material{end+1} = stationIB.sm{kdp};
                case {'flare','hourglass'}
                    % flare and hourglass are equivalent on the area inboard edge
                    % start and end of first area is current DP
                    % start and end of next area is current and next DP
                    SkinAreas(kStation).startIB(end+(1:2)) = [kdp kdp];
                    SkinAreas(kStation).endIB(end+(1:2))   = [kdp kdp+1];
                    SkinAreas(kStation).Material{end+1} = stationIB.dpmaterial{kdp};  %jcb:  append (end+1) dp material
                    SkinAreas(kStation).Material{end+1} = stationIB.sm{kdp};  %jcb:  append (end+1) segment material
            end
        end

        for kdp = 1:numel(stationOB.dp)-1
            switch stationOB.dptype{kdp}
                case {'single','flare'}
                    % single and flare are equivalent on the area outboard edge
                    % start and end are current and next DP
                    SkinAreas(kStation).startOB(end+1) = kdp;
                    SkinAreas(kStation).endOB(end+1)   = kdp+1;
                case {'double','hourglass'}
                    % double and hourglass are equivalent on the area outboard edge
                    % start and end of first area is current DP
                    % start and end of next area is current and next DP
                    SkinAreas(kStation).startOB(end+(1:2)) = [kdp kdp];
                    SkinAreas(kStation).endOB(end+(1:2))   = [kdp kdp+1];
            end
        end
    end
    %tcl: Determine which composite materials are used in the model
    compsInModel = {};
    %tcl: search shear web materials
    for k = 1:TotalShearwebs
        compsInModel = [compsInModel; {app.shearweb(k).Material}];  %#ok
    end
    %tcl: search skin materials
    for k = 1:numel(SkinAreas)
        compsInModel = [compsInModel; transpose(SkinAreas(k).Material)];  %#ok
    end
    compsInModel = unique(compsInModel);

    % load the material database and create searchable list
    if ~isempty(app.settings.job_name)
        % job name exists, use the local material database
        app.matdb_path = fullfile(app.settings.job_path,'MatDBsi.txt');
    else
        % if no job name, use the master material database
        % jcb: we shouldn't arrive here because a file save is required first
        app.matdb_path = fullfile(app.numadpath,'MatDBsi.txt');
    end
    app.matdb = readMatDB(app.matdb_path);
    for k=1:numel(app.matdb)
        app.matlist{k} = app.matdb(k).name;
        mattype{k} = app.matdb(k).type;  %#ok
    end
    app.isotropic =  strcmp('isotropic',mattype);
    app.orthotropic = strcmp('orthotropic',mattype);
    app.composite = strcmp('composite',mattype);

    % Determine which isotropic and orthotropic materials are used in the model
    isoorthoInModel = {};
    for kcomp = 1:numel(compsInModel)
        n = strcmp(compsInModel(kcomp),app.matlist);
        if ~any(n)
            errordlg(sprintf('Material "%s" not found in database.',compsInModel{kcomp}),'Error');
            error('Material "%s" not found in database.',compsInModel{kcomp});
        end
        mat = app.matdb(n);
        layerNames = cell(numel(mat.layer),1);
        for klay = 1:numel(mat.layer)
            layerNames(klay) = {mat.layer(klay).layerName};
        end
        isoorthoInModel = unique([isoorthoInModel; layerNames]);
    end
