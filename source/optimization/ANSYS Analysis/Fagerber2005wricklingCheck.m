%% Wrinkling Analysis Post-Processing Script
function [limitingElementData]=Fagerber2005wricklingCheck(app,SkinAreas,compsInModel,coreMatName)
    %limitingElementData - [ansysSecNumber elno lf phicr]
    TotalStations = numel(app.station);
    TotalShearwebs = numel(app.shearweb);

    %%%%%%%%%%%%%%%%%   Main loop #1: loop around aero shell.   %%%%%%%%%%%%%%%%%
    LF=[];
    for kStation = 1:TotalStations-1 %Loop along span
        for kArea = 1:numel(SkinAreas(kStation).startIB)     
            %See if the section contatins Balsa/core material name (i.e find
            %the sandwhich panels)
            n = strcmp(SkinAreas(kStation).Material{kArea},app.matlist);
            mat = app.matdb(n);
            if contains([mat.layer.layerName],coreMatName)  %%%%%%%%%%%%%%%%%This logic will need to be made more general
                ansysSecNumber = find(strcmp(SkinAreas(kStation).Material(kArea),compsInModel)==1);
                file=strcat('section-',int2str(ansysSecNumber),'-faceAvgStresses.txt');
                avgFaceStress=txt2mat(file);  
                delete(file)
                LF=getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber,avgFaceStress,app,mat,coreMatName);
            end
        end
    end

     
    %%%%%%%%%%%%%%%%%   Main loop #2: loop along web.   %%%%%%%%%%%%%%%%%
    for kShearweb = 1:TotalShearwebs
        n = strcmp(app.shearweb(kShearweb).Material,app.matlist);
        mat = app.matdb(n);
        if contains([mat.layer.layerName],coreMatName)  %%%%%%%%%%%%%%%%%This logic will need to be made more general
            ansysSecNumber = find(strcmp({app.shearweb(kShearweb).Material},compsInModel)==1);
            file=strcat('section-',int2str(ansysSecNumber+1000),'-faceAvgStresses.txt');
            avgFaceStress=txt2mat(file);
            delete(file)
            LF=getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber+1000,avgFaceStress,app,mat,coreMatName);
        end
    end

    [minLF,index]=min(LF(:,3));
    
    limitingElementData=LF(index,:);
    fprintf('\n\n The minimum wrinkling LF is: %f, wrinkle angle: %.2f°' ,minLF, LF(index,4))
    fprintf('\n and occurs in section number %i, element number %i\n, ',LF(index,1),LF(index,2))

    % [maxLF,index]=max(LF(:,3));
    % fprintf('\n\n The maximum LF is: %f, wrinkle  angle: %.2f°' ,maxLF, LF(index,4))
    % fprintf('\n and occurs in section number %i, element number %i, ',LF(index,1),LF(index,2))

end