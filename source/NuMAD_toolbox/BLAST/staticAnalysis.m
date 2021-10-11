function [displ,staticAnalysisSuccessful]=staticAnalysis(model,displ,Omega,OmegaStart,smOmega,aeroLoadsFlag,analysisType)
%staticAnalysis Performs nonlinear static elasticity analysis
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [displ,staticAnalysisSuccessful]=staticAnalysis(outfile,model,displ,
%                                    Omega,OmegaStart,smOmega,aeroLoadsFlag,
%                                    analysisType)
%                    
%   This function performs nonlinear static analysis through a direct
%   iteration procedure. 
%
%      input:
%      outfile          = output filename
%      model            = model struct
%      displ            = displacement vector (intial guess)
%      Omega            = rotor speed (Hz)
%      OmegaStart       = start rotor speed (for load stepping)
%      smOmega          = mode natural frequency for calculation of
%                         unsteady aero loads
%      aeroLoadsFlag    = flag to include or exclude aerodynamic loading
%                         0-exclude, 1-include
%      analysisType     = char for stability analysis type
%      
%      output:
%      displ                     = converged displacement output
%      staticAnalysisSuccessful  = boolean describing if static analysis was
%                                  successful (converged)


tic

%initialize variables
numEl = model.numElements;
x = model.mesh.x;
y = model.mesh.y;
z = model.mesh.z;
conn = model.mesh.conn;
numNodes = length(x);
elementOrder = model.elementOrder;
BC = model.BC;

numNodesPerEl = elementOrder + 1;
numDOFPerNode = 6;
totalNumDOF = numNodes * numDOFPerNode;

Kg = zeros(totalNumDOF,totalNumDOF);
Fg = zeros(totalNumDOF,1);

%define load steps for nonlinear analysis
[loadSteps] = defineLoadSteps(Omega,OmegaStart,3600.0);
numLoadSteps=length(loadSteps);

%define max iterations before terminating nonlinear analysis
MAXIT=1000;

%load step loop
for m=1:numLoadSteps
staticAnalysisSuccessful = false;
uNorm = 1.0;

iterationCount = 0;
while(uNorm > 1.0e-6 && iterationCount < MAXIT)
    Kg = zeros(totalNumDOF,totalNumDOF);
    Fg = zeros(totalNumDOF,1);

    for i=1:numEl
        
        %update solution from last step (with convergence acceleration methods)
        if(iterationCount >=1)
            gamm = 0.5;
            dispEval = dispOld.*gamm + displ.*(1-gamm);
        else
            dispEval = displ;
        end
        
        %Calculate Ke and Fe for element i
        
        %determine element length
        elxLoc = [0.0 model.mesh.elLength(i)];  
        
        index = 1;
        for j=1:numNodesPerEl
            elx(j) = x(conn(i,j));
            ely(j) = y(conn(i,j));
            elz(j) = z(conn(i,j));
            for k=1:numDOFPerNode
                eldisp(index) = dispEval((conn(i,j)-1)*numDOFPerNode + k);
                index = index + 1;
            end
        end
        
        %calculate element stiffness matrix and load vector (includes spin
        %softening)
        [Ke,~,~,Fe] = calculateEulerBeamElement(elementOrder,elx,ely,elz,elxLoc,model.hubRadius,eldisp,model.sectionPropsArray{i},model.sweepAngle(i),model.coneAngle(i),model.aeroSweepAngle(i),Omega,smOmega,model.airDensity,aeroLoadsFlag,false,analysisType);
      
        %Assemble element i into global system
        [Kg,Fg] = assembly(Ke,Fe,conn(i,:),numNodesPerEl,numDOFPerNode,Kg,Fg); 

    end
    
    %scale load vector based off of load steps
    Fg=loadSteps(m).*Fg;
    %----------------------------------------------------------------------
    %APPLY BOUNDARY CONDITIONS
    [Kg,Fg] = applyBC(Kg,Fg,model.BC,dispEval,0,numDOFPerNode);

    %set solution at previous iteration
    dispOld = displ;

    %solve for current iteration solution
    displ = Kg\Fg;
    
    %calculate norm of difference between current and previous iteration
    %solution
    uNorm = calculateNorm(dispOld,displ);
    iterationCount = iterationCount +1;
end

    if(iterationCount>=MAXIT)
        disp('Max iterations exceeded');
        break; %breaks out of load step loop
    else
        staticAnalysisSuccessful=true;
    end
end


    
    t_static = toc;
    
end


function [uNorm] = calculateNorm(uPrev,u)
% This function calculates the relative norm between two vectors
%
%    input:
%    uPrev   = displacement vector from previous step
%    u       = displacement vector from current step
%
%    output:
%    uNomrm  = calculated relative norm of u - uPrev

len = length(u);
    numerator = 0;
    denom = 0;
    for i=1:len
        numerator = numerator + (u(i)-uPrev(i))^2;
        denom = denom + u(i)^2;
    end
    
    uNorm = sqrt(numerator/denom);
end
% 
% function dispPlotter(q,dofPerNode)
% 
% len = length(q);
% len2 = len/dofPerNode;
% qmat = zeros(dofPerNode,len2);
% 
% for i=1:len2
%     index = dofPerNode*(i-1) +1;
%     qmat(:,i) = q(index:index+(dofPerNode-1));
% end
% 
% for i=1:dofPerNode
% figure(i);
% plot(qmat(i,:));
% end
% 
% end

function [loadSteps] = defineLoadSteps(maxOmega,OmegaStart,minLoadStep)
% This function defines load steps
%
%    input:
%    maxOmega    = maximum omega value
%    OmegaStart  = starting omega value
%    minLoadStep = minimum desires load step
%
%    output:
%    loadSteps   = calculated load steps

temp = [OmegaStart^2+minLoadStep:minLoadStep:maxOmega^2];

if(isempty(temp))
    temp=maxOmega^2;
end

len=length(temp);

if(temp(len)~=maxOmega^2)
    temp = [temp, maxOmega^2];
end
if(temp(1)==0.0 && maxOmega==0.0)
    loadSteps = 1.0;
else
    loadSteps = temp./(maxOmega^2);
end
end
