function angleChange=getAirfoilNormalsAngleChange(unitNormals)
%Find the angle changes between adjacent unit vectors
    angleChange=zeros(length(unitNormals),1);
    for iVector =1:length(unitNormals)-1
        currentVector=unitNormals(iVector,:);
        nextVector=unitNormals(iVector+1,:);
        angleChange(iVector)=acosd(currentVector*nextVector');
    end

    %angle change between last point and first point
    currentVector=unitNormals(iVector+1,:);
    nextVector=unitNormals(1,:);
    angleChange(iVector+1)=acosd(currentVector*nextVector');
end