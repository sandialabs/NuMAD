function [out] = loadFASTOutDataGageRot(fileName,gageCoordinateRotation)
    %% Get the data from a particular fast file, with the forces and moments rotated
    %% to account for structural twist, based on gageCoordinateRotation.
    out=loadFASTOutData(fileName);  % load a file from the list
    
    % correct for spanwise gage channels which are in the local blade
    % coordinate system oriented along the principal axis
    for gg = 1:length(gageCoordinateRotation)
        % determine the rotation matrix
        theta = gageCoordinateRotation(gg);
        rotMat = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
        
        % force transformation
        fxData = out.data(:,strcmp(['Spn' num2str(gg) 'FLxb1'],out.list));
        fyData = out.data(:,strcmp(['Spn' num2str(gg) 'FLyb1'],out.list));     
        Ftransform = rotMat*[fxData'; fyData'];
        out.data(:,strcmp(['Spn' num2str(gg) 'FLxb1'],out.list)) = Ftransform(1,:)';
        out.data(:,strcmp(['Spn' num2str(gg) 'FLyb1'],out.list)) = Ftransform(2,:)';
        out.tbl.(['Spn' num2str(gg) 'FLxb1']) = Ftransform(1,:)';
        out.tbl.(['Spn' num2str(gg) 'FLyb1']) = Ftransform(2,:)';
        
        % moment transformation
        mxData = out.data(:,strcmp(['Spn' num2str(gg) 'MLxb1'],out.list));
        myData = out.data(:,strcmp(['Spn' num2str(gg) 'MLyb1'],out.list));
        Mtransform = rotMat*[mxData'; myData'];
        out.data(:,strcmp(['Spn' num2str(gg) 'MLxb1'],out.list)) = Mtransform(1,:)';
        out.data(:,strcmp(['Spn' num2str(gg) 'MLyb1'],out.list)) = Mtransform(2,:)';
        out.tbl.(['Spn' num2str(gg) 'MLxb1']) = Mtransform(1,:)';
        out.tbl.(['Spn' num2str(gg) 'MLyb1']) = Mtransform(2,:)';
    end
end

