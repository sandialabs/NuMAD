function unitNormals=getAirfoilNormals(coordinates)
    %Method finds which airfoil is flatback. If points are placed
    %in flatback region, they are removed for good resampling
    %results. Currently this removal only works
    %for one point located on TE region. Method also gets the TE tpye for round sections. 
    
    %coordinates m by 2 matrix where m is the number of points about
    %airfoil. Coordiantes are 2D.
    
    nPoints=length(coordinates);
    unitNormals=zeros(nPoints-1,2);
    for iPoint=1:nPoints-1
       x1=coordinates(iPoint,1);
       y1=coordinates(iPoint,2);
       x2=coordinates(iPoint+1,1);
       y2=coordinates(iPoint+1,2);
       currentPoint = [x1, y1]; nextPoint = [x2, y2]; 
       xCoords=[x1;x2]; yCoords=[y1;y2];
       r=nextPoint-currentPoint; %Postion vector from currentPoint to nextPoint
       if abs(r(1)) + abs(r(2)) ~= 0 %Skip if points are coincedint
           unitNorm=null(r)';
           crossProduct=cross([r 0],[unitNorm 0]);
           if crossProduct(3)<0
               unitNorm=-unitNorm;
           end
           unitNormals(iPoint,:)=unitNorm;
       else
           unitNormals(iPoint,:)=[NaN NaN];
       end
    end
end
%             for iPoint=1:nPoints
%                 text(coordinates(iPoint,1),coordinates(iPoint,2),num2str(iPoint),'Color','b')
%             end
