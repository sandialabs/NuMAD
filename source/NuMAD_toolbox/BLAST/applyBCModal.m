function [K,dofVector] = applyBCModal(K,BC,numDofPerNode)
%applyBCModal Applies boundary conditions to system for modal analysis
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [K,dofVector] = applyBCModal(K,BC,numDofPerNode)
%                    
%   This function applies boundary conditions to a system matrix for modal
%   analysis
%
%      input:
%      K             = assembled global system matrix
%      BC            = struct of boundary condition information
%      numDofPerNode = number of degrees of freedom per node
 
%      output:
%      K             = global system matrix with boundary conditions
%      dofVector     = reduced DOF vector after imposing BCs




[numEq,dum]=size(K);

[dofVector] = calculateReducedDOFVector(numEq,numDofPerNode,BC);

%APPLY BCs FOR PRIMARY VARIABLE
if(BC.numpBC > 0)
    pBC = BC.pBC;
    [numpBC,dum] = size(pBC);
    
    for i=1:numpBC
        nodeNumber = pBC(i,1);
        dofNumber = pBC(i,2);
        specVal = pBC(i,3);
        
        eqNumber = (nodeNumber-1)*numDofPerNode + dofNumber-(i-1);
        
        for j=eqNumber+1:numEq
            K(j-1,:) = K(j,:);
            K(:,j-1) = K(:,j);
        end
    end
end

%resize matrix
K=K(1:numEq-numpBC,1:numEq-numpBC);
end


function [dofVector] = calculateReducedDOFVector(totalDOFs,numDofPerNode,BC)

    numNodes = totalDOFs/numDofPerNode;
    index = 1;
    for i=1:numNodes
       for j=1:numDofPerNode
           [constrained] = checkIfConstrained(i,j,BC);
           if(constrained == false)
            dofVector(index,1) = (i-1)*numDofPerNode + j;
            dofVector(index,2) = i;
            dofVector(index,3) = j;
            index = index + 1;
           end
       end
    end

end

function [constrained] = checkIfConstrained(nodeNum,localDofNum,BC)
    constrained = false;
    for i = 1:BC.numpBC
       if ((BC.pBC(i,1) == nodeNum) && (BC.pBC(i,2) == localDofNum))
           constrained = true;
           break;
       end
    end
end