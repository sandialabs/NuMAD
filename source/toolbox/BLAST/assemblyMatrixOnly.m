function [Kg] = assemblyMatrixOnly(Ke,conn,numNodesPerEl,numDOFPerNode,Kg) 
%assemblyMatrixOnly Assembles element matrices into global sys of equations
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [Kg] = assemblyMatrixOnly(Ke,conn,numNodesPerEl,numDOFPerNode,Kg) 
%                    
%   This function assembles the element matrix into the
%   global system of equations
%
%      input:
%      Ke             = element matrix
%      conn           = element connectivity
%      numNodesPerEl  = number of nodes per element
%      numDofPerNode  = number of degrees of freedom per node
%      Kg             = global system matrix
 
%      output:
%      Kg             = global system matrix with assembled element

count = 1;
for i=1:numNodesPerEl
   for j=1:numDOFPerNode
       dofList(count) = (conn(i)-1)*numDOFPerNode + j;
       count = count +1;
   end
end

numDOFPerEl = length(dofList);
%Assemble element i into global system
        for j=1:numDOFPerEl
            J = dofList(j);
            for m=1:numDOFPerEl
                M = dofList(m);
                Kg(J,M) = Kg(J,M) + Ke(j,m);
            end
        end
end
        