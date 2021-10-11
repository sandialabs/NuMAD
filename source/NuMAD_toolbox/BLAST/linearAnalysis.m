function  [freq,damp,phase1,phase2,imagCompSign,convergedCheck] = linearAnalysis(model,displ,Omega,smOmega,aeroLoadsFlag,analysisType)
%linearAnalysis Performs modal structural dynamics analysis
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [displ,staticAnalysisSuccessful]=staticAnalysis(outfile,model,displ,
%                                    Omega,OmegaStart,smOmega,aeroLoadsFlag,
%                                    analysisType)
%                    
%   This function performs linear modal structural dynamics analysis.
%
%      input:
%      outfile          = output filename
%      model            = model struct
%      displ            = displacement vector
%      Omega            = rotor speed (Hz)
%      OmegaStart       = start rotor speed (for load stepping)
%      smOmega          = mode natural frequency for calculation of
%                         unsteady aero loads
%      aeroLoadsFlag    = flag to include or exclude aerodynamic loading
%                         0-exclude, 1-include
%      plotFlag         = flag to plot mode shapes 0-off, 1-on
%      analysisType     = char for stability analysis type
%      
%      output:
%      freq             = frequencies of calculated modes
%      damp             = damping ratio of calculate modes
%      modesAmp         = amplitudes of calculated modes
%      modesPhase       = phase of calculated modes
%      imagCompSign     = sign of imaginary component on eigenvalue

tic

% intialization of variables
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
Mg = zeros(totalNumDOF,totalNumDOF);
Cg = zeros(totalNumDOF,totalNumDOF);

modesAmp = 0;
modesPhase=0;
imagCompSign=0;

%begin element calculation loop
    for i=1:numEl
        %define element length
        elxLoc = [0.0 model.mesh.elLength(i)];
        
        %get element coordinates and displacement values
        index = 1;
        for j=1:numNodesPerEl
            elx(j) = x(conn(i,j));
            ely(j) = y(conn(i,j));
            elz(j) = z(conn(i,j));
            for k=1:numDOFPerNode
                eldisp(index) = displ((conn(i,j)-1)*numDOFPerNode + k);
                index = index + 1;
            end
        end
        %calculate Ke, Me, Ce for element i
        [Ke,Me,Ce,~] = calculateEulerBeamElement(elementOrder,elx,ely,elz,elxLoc,model.hubRadius,eldisp,model.sectionPropsArray{i},model.sweepAngle(i),model.coneAngle(i),model.aeroSweepAngle(i),Omega,smOmega,model.airDensity,aeroLoadsFlag,true,analysisType);

        %Assemble element i into global system
        [Kg] = assemblyMatrixOnly(Ke,conn(i,:),numNodesPerEl,numDOFPerNode,Kg); 
        [Mg] = assemblyMatrixOnly(Me,conn(i,:),numNodesPerEl,numDOFPerNode,Mg); 
        [Cg] = assemblyMatrixOnly(Ce,conn(i,:),numNodesPerEl,numDOFPerNode,Cg); 
    end
    
    %----------------------------------------------------------------------
    %APPLY BOUNDARY CONDITIONS
    [Kg,dofVector] = applyBCModal(Kg,model.BC,numDOFPerNode);
    [Mg,dofVector] = applyBCModal(Mg,model.BC,numDOFPerNode);
    [Cg,dofVector] = applyBCModal(Cg,model.BC,numDOFPerNode);
    
    %set solver flag
    if(Omega==0.0)
        solveFlag = 1;
    else
        solveFlag = 2;
    end
    
    if(strcmpi(analysisType,'D'))
        solveFlag = 2;
        Cg=Cg.*0;
    end
    
    %perform eigensolve on system
    [eigVec,eigVal] = eigSolve(Mg,Cg,Kg,solveFlag);
    
    %process eigenvalues and eigenvectors for frequency, damping, and mode
    %shapes
    
    %initialize post-process variables
    [~,len] = size(eigVal);
    eigValRe = zeros(len,1);
    eigValIm = eigValRe;
    freq = eigValRe;
    damp = eigValRe;
%     phase1 = eigValRe;
%     phase2 = eigValRe;
    imagCompSign = eigValRe;
    
    for i=1:len
        eigValRe(i) = real(eigVal(i,i));
        eigValIm(i) = imag(eigVal(i,i));
        [freq(i),damp(i),phase1(:,:,i),phase2(:,:,i)] = extractFreqDamp(eigVal(i,i),eigVec(:,i),numDOFPerNode,dofVector,analysisType);
        imagCompSign(i) = sign(imag(eigVal(i,i)));
    end
    
    %check frequencies for convergence
    convFreqTol = 0.01;
    [convergedCheck] = checkFrequenciesForConvergence(freq,smOmega,convFreqTol);
    
    %sort modes
    [freq,damp,phase1,phase2,imagCompSign,convergedCheck] = ...
       sortModes(freq,damp,phase1,phase2,imagCompSign,convergedCheck);
     
end



function [uNorm] = calculateNorm(uPrev,u)

%This function calculate the norm of the difference between two vectors.

%    input:
%    uPrev = vector of displacements from last time step
%    u     = vector of displacements from current timestep
%
%    output:
%    uNorm = relative norm of the difference between uPrev and u
%
len = length(u);
    numerator = 0;
    denom = 0;
    for i=1:len
        numerator = numerator + (u(i)-uPrev(i))^2;
        denom = denom + u(i)^2;
    end
    
    uNorm = sqrt(numerator/denom);
end

function [convergedCheck] = checkFrequenciesForConvergence(freq,guessFreq,tol)

%This function checks frequencies for convergence

%    input:
%    freq      = calculated frequencies
%    guessFreq = guess frequency from p-k iteration procedure
%    tol       = tolerance used to define convergence
%
%    output:
%    convergedCheck = list of booleans denoting if entries in freq are
%                     converged
%
len = length(freq);
guessFreq = ones(len,1)*guessFreq;
convergedCheck(1:len) = false;

del = abs(freq-guessFreq)./guessFreq;

for i=1:len
   if(del(i)<tol)
       convergedCheck(i) = true;
   end
end

end
