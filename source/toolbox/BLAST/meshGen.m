function [mesh] = meshGen(numElements,secLocations,hubRadius,preSwp,preCrv,edgEAOff,flpEAOff,elementOrder,mesh)
%meshGen  Generates finite element  beam mesh of HAWT blade for BLAST
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [mesh] = meshGen(numElements,secLocations,hubRadius,preSwp,edgEAOff,
%                    elementOrder,mesh)
%   This function accepts mesh discretization information and geometry 
%   obtained from the blade file to construct the finite element mesh
%
%      input:
%      numElements  = number of elements in generated mesh
%      secLocations = radial location of blade section locations
%      hubRadius    = hub radius of turbine rotor (offset of blade root 
%                     from axis of rotation)
%      preSwp       = preSweep parameter from blade file at section
%                     locations to characterize blade sweep
%      edgEAOff     = edgewise offset of elastic axis from reference axis
%                     at section locations 
%      elementOrder = order of elements to be constructed, alters number of
%                     nodes per element. This option is written for future 
%                     flexibility. Only linear elements have been 
%                     implemented.
%      mesh         = mesh struct input containing edgewise nodal
%                     coordinate read in from earlier stage in the code
%
%      output:
%      mesh         = mesh output struct containing x, y, z nodal
%                     coordinates, element length, and element connectivity
%

    x = secLocations+hubRadius;
    y = preSwp + edgEAOff;
    z = preCrv + flpEAOff;
    
    %determine connectivity
    if(elementOrder == 1)
        conn(1,:)=[1,2];
    elseif(elementOrder == 2)
        conn(1,:)=[1,2,3];
    end
    for i=2:numElements
        conn(i,1) = conn(i-1,1+elementOrder);
       for j=1:elementOrder
            conn(i,1+j) = conn(i,j)+1;
       end
    end
    
    mesh.x=x;
    mesh.y=y;
    mesh.z=z;

    for i=1:numElements
        delX=mesh.x(i+1)-mesh.x(i);
        delY=mesh.y(i+1)-mesh.y(i);
        delZ=mesh.z(i+1)-mesh.z(i);
        mesh.elLength(i)=sqrt(delX^2 + delY^2 + delZ^2);
    end
    
    mesh.conn=conn;
    
end