function [Kg,Fg] = assembly(Ke,Fe,conn,numNodesPerEl,numDOFPerNode,Kg,Fg) 
%assembly Assembles element matrices into global system of equations
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [Kg,Fg] = assembly(Ke,Fe,conn,numNodesPerEl,numDOFPerNode,Kg,Fg) 
%                    
%   This function assembles the element matrix and load vector into the
%   global system of equations
%
%      input:
%      Ke             = element matrix
%      Fe             = element vector
%      conn           = element connectivity
%      numNodesPerEl  = number of nodes per element
%      numDofPerNode  = number of degrees of freedom per node
%      Kg             = global system matrix
%      Fg             = global load vector
 
%      output:
%      Kg             = global system matrix with assembled element
%      Fg             = global load vector with assembled element

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
            Fg(J) = Fg(J) + Fe(j);
            for m=1:numDOFPerEl
                M = dofList(m);
                Kg(J,M) = Kg(J,M) + Ke(j,m);
            end
        end
end
        