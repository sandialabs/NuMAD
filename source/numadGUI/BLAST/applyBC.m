function [Kg,Fg] = applyBC(Kg,Fg,BC,u,iterationType,numDofPerNode)
%applyBC Applies boundary conditions to system for static analysis
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [Kg,Fg] = applyBC(Kg,Fg,BC,u,iterationType,numDofPerNode)
%                    
%   This function applies boundary conditions to the stiffness matrix and
%   load vector for a static analysis.
%
%      input:
%      Kg            = assembled global stiffness matrix
%      Fg            = assembled global load vector
%      BC            = struct of boundary condition information
%      u             = global displacement vector
%      iterationType = for nonlinear analysis, not used in BLAST
%      numDofPerNode = number of degrees of freedom per node
 
%      output:
%      Kg            = global stiffness matrix with boundary conditions
%      Fg            = global load vector with boundary condition


[numEq,dum]=size(Kg);

%APPLY BCs FOR PRIMARY VARIABLE
if(BC.numpBC > 0)
    pBC = BC.pBC;
    [numpBC,dum] = size(pBC);
    
    for i=1:numpBC
        nodeNumber = pBC(i,1);
        dofNumber = pBC(i,2);
        specVal = pBC(i,3);
        
        eqNumber = (nodeNumber-1)*numDofPerNode + dofNumber;
        
        for j=1:numEq
            Kg(eqNumber,j) = 0.0;
            Fg(j) = Fg(j) - Kg(j,eqNumber)*specVal;
            Kg(j,eqNumber) = 0.0;
        end
        Fg(eqNumber) = specVal;
        Kg(eqNumber,eqNumber) = 1.0;
    end
end

%APPLY BCs FOR SECONDARY VARIABLE
if(BC.numsBC > 0)
    sBC = BC.sBC;
    [numsBC,dum] = size(sBC);
    
    for i=1:numsBC
       nodeNumber = sBC(i,1);
       dofNumber = sBC(i,2);
       specVal =  sBC(i,3);
       
       eqNumber = (nodeNumber-1)*numDofPerNode + dofNumber;
       
       Fg(eqNumber) = Fg(eqNumber) + specVal;
        
    end
end

%APPLY BCs FOR MIXED BCs

if(BC.nummBC > 0)
    mBC = BC.mBC;
    [nummBC,dum] = size(mBC);
    
    for i=1:nummBC
       nodeNumber = mBC(i,1);
       dofNumber  = mBC(i,2);
       beta0 =  mBC(i,3);
       beta1 =  mBC(i,4);
       uRef   =  mBC(i,5);
              
       eqNumber = (nodeNumber-1)*numDofPerNode + dofNumber;
       uVal   = u(eqNumber);
       if(iterationType<=1)
            Fg(eqNumber) = Fg(eqNumber) + (beta0 + beta1*uVal)*uRef;
            Kg(eqNumber,eqNumber) = Kg(eqNumber,eqNumber) + (beta0 + beta1*uVal);
       else
           Fg(eqNumber) = Fg(eqNumber) - (beta0 + beta1*uVal)*uVal + (beta0+beta1*uVal)*uRef;
           Kg(eqNumber,eqNumber) = Kg(eqNumber,eqNumber) + (beta0 + 2*beta1*uVal)-beta1*uRef;
       end
    end
end
end