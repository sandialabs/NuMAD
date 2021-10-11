clear all
switch 1
    case 1
        pn = 'C:\Users\jcberg\Documents\NuMAD\NuMAD_models\CX100_v1.0 - Copy\';
        fn = 'CX100_v1.0.nmd';
        maxLayers = 18;
    case 2
        pn = 'C:\Users\jcberg\Documents\NuMAD\NuMAD_models\SNLWindPACTv1.3\';
        fn = 'WindPact15MWv1.3.nmd';
        maxLayers = 6;
    case 3
        pn = 'C:\Users\jcberg\Documents\NuMAD\NuMAD_models\SNL61p5m_NuMAD_121212\';
        fn = 'SNL61p5m.nmd';
        maxLayers = 120;
end
input_file = fullfile(pn,fn);
app.userpath = getenv('APPDATA');
numadfolder = 'NuMAD';
app.userpath = fullfile(app.userpath,numadfolder);
app.settings = readNuMADsettings(fullfile(app.userpath,'settings.txt'));
app.settings.job_name = fn;
app.settings.job_path = pn;
[app.station app.shearweb app.active app.ansys app.BladeRotation blade app.plot3d app.flutterInput] = readNuMADinput(input_file);
app.n_panels = 200;
app.station = resampleAirfoilDB(app.station,app.n_panels,'cosine');
[app.afdb app.aflist] = readAirfoilDB(fullfile(pn,'airfoils'));  % use local database
app.afdb = resampleAirfoilDB(app.afdb,app.n_panels,'cosine');
app.matdb = readMatDB(fullfile(pn,'MatDBsi.txt'));  % use local database

%%
[KeyPoints, SkinAreas, isoorthoInModel] = write_shell7(app,blade);
matlist = {app.matdb.name};
colors = cool(numel(isoorthoInModel));

%% getStacks()
nStations = numel(SkinAreas);
nStacks = 0;
for k=1:nStations
    nStacks = max(nStacks,numel(SkinAreas(k).Material));
end
nUniqLayers = 0;
for k=1:numel(app.matdb)
    if isempty(app.matdb(k).layer)
        continue
    end
    nUniqLayers = max(nUniqLayers,numel(app.matdb(k).layer));
end
StackArray = zeros(nStations,nUniqLayers,nStacks);
for kStation = 1:nStations
    for kStack = 1:numel(SkinAreas(kStation).Material)
        n = find(strcmp(SkinAreas(kStation).Material{kStack},matlist)==1);
        mat = app.matdb(n);
        for kLayer = 1:numel(mat.layer)
            n = find(strcmp(mat.layer(kLayer).layerName,matlist)==1);
            StackArray(kStation,kLayer,kStack) = n;
        end
    end
end


%%
figure(100); clf; axis off;
for kColor = 1:length(colors)
    x = [1 10 10 1];
    y = [0 0  1  1] + kColor;
    patch(x,y-0.5,colors(kColor,:));
    text(x(1)+1,y(1),isoorthoInModel(kColor),'interpreter','none');
end
%%
figure(1); clf;
for layerToPlot = 1:6 %1:maxLayers
subplot(2,3,layerToPlot)
title(sprintf('%s, layer=%d',fn,layerToPlot),'interpreter','none')
for kStation = 1:numel(SkinAreas)
    for kArea = 1:numel(SkinAreas(kStation).Material)
        n = find(strcmp(SkinAreas(kStation).Material{kArea},matlist)==1);
        mat = app.matdb(n);
        matLayerIndices = [];
        for kLayer = 1:numel(mat.layer)
%             matLayerIndices = [matLayerIndices, repmat(kLayer,1,mat.layer(kLayer).quantity)];
            matLayerIndices = [matLayerIndices, kLayer];
        end
        if layerToPlot <= numel(matLayerIndices)
            kLayer = matLayerIndices(layerToPlot);
            n = find(strcmp(mat.layer(kLayer).layerName,isoorthoInModel)==1);
            kcolor = rem(n-1,numel(colors))+1;
            patchColor = colors(kcolor,:);
        else
            patchColor = [1 1 1]; %white
        end
        ib.x = kStation;
        ob.x = kStation+1;
        ib.y1 = SkinAreas(kStation).startIB(kArea);
        ob.y1 = SkinAreas(kStation).startOB(kArea);
        ob.y2 = SkinAreas(kStation).endOB(kArea);
        ib.y2 = SkinAreas(kStation).endIB(kArea);
        x = [ib.x  ob.x  ob.x  ib.x ];
        y = [ib.y1 ob.y1 ob.y2 ib.y2];
        hg_area(kStation,kArea) = patch(x,y,patchColor);
    end
end
% pause(0.5)
% pause
end