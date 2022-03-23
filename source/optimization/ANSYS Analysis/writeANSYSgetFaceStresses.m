%% Panel Stresses Analysis Script
function [app,SkinAreas,compsInModel]=writeANSYSgetFaceStresses(blade,fid,coreMatName)
  
    [isoorthoInModel,compsInModel,SkinAreas,app] = getMatrialLayerInfoWithOutGUI(blade);
    %fid=fopen('getFaceStresses.mac','w+');
    TotalStations=numel(blade.ispan);
    for kStation = 1:TotalStations-1 %Loop along span
        %kPanel=find(~cellfun('isempty',strfind([SkinAreas(kStation).Material],'PANEL'))); %Array that stores the kArea index that contains 'PANEL' in the name
        %for i=1:numel(kPanel)
        for kArea = 1:numel(SkinAreas(kStation).startIB)

            %See if the section contatins Balsa/core material name (i.e find
            %the sandwhich panels)
            n = strcmp(SkinAreas(kStation).Material{kArea},app.matlist);
            mat = app.matdb(n);
            if contains([mat.layer.layerName],coreMatName)  %%%%%%%%%%%%%%%%%This logic will need to be made more general
                ansysSecNumber = find(strcmp(SkinAreas(kStation).Material(kArea),compsInModel)==1);

                writeANSYSinputFile(fid,mat,ansysSecNumber,coreMatName)

            end
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'!*************** WEB ***************\n');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    %fclose(fid);
    %Web
    TotalShearwebs = numel(app.shearweb);
    for kShearweb = 1:TotalShearwebs
        n = strcmp(app.shearweb(kShearweb).Material,app.matlist);
        mat = app.matdb(n);
        if contains([mat.layer.layerName],coreMatName)  %%%%%%%%%%%%%%%%%This logic will need to be made more general
            ansysSecNumber = find(strcmp({app.shearweb(kShearweb).Material},compsInModel)==1);
            ansysSecNumber=ansysSecNumber+1000;
            writeANSYSinputFile(fid,mat,ansysSecNumber,coreMatName)

        end
    end
    fprintf(fid,'FINISH\n');
    fprintf(fid,'allsel\n');

end