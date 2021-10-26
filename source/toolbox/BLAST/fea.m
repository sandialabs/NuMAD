function fea(inStruct,Omega,smOmega,aeroLoadsFlag,analysisType)
%fea Performs flutter analysis for single rotor speed and freq guess
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   fea(infile,outfile,Omega,smOmega,aeroLoadsFlag,analysisType)
%
%   This function reads input, generates a beam element mesh, performs
%   static "spin up" analysis, and modal analysis of system at a given
%   rotor speed and frequency guess.
%
%      input:
%      inStruct      = struct with input data
%      Omega         = rotor speed (RPM)
%      smOmega       = flutter frequency guess
%      aeroLoadsFlag = flag to include or exclude aerodynamic loads
%                      (1-include, 0-exclude)
%      analysisType  = char denoting stability analysis type (F-rotating
%                      flutter analysis, P-parked flutter analysis, 
%                      D-rotating divergence analysis)
%      
%      output:
%      none
%      no explicit output is returned from this function but output data is
%      written to outfile

%Set analysis type
if(strcmp(analysisType,'F'))
    disp('Running in rotating flutter mode');
elseif(strcmp(analysisType,'P'))
    disp('Running in parked flutter mode');
elseif(strcmp(analysisType,'D'))
    disp('Running rotating divergence mode');
else
    error('Analysis type not recognized. Exiting.');
end

tic

%open input and output files, read in input
fid=fopen(inStruct.outFile,'wt');

fstfn=inStruct.fstFile;
bladefn = inStruct.bladeFile;
adfn= inStruct.aeroFile;
pitchAxisLoc = inStruct.pitchAxisDomain;
pitchAxisMat = inStruct.pitchAxisVal;
LCSvalsMat = inStruct.LCS;

model = readInput(fstfn,bladefn,adfn,pitchAxisLoc,pitchAxisMat,LCSvalsMat);
t_input = toc;
disp('Elapsed time for input(s):');
disp(t_input);

% %------------------------------

tic
%-------- MESH GENERATION -----
model.mesh = meshGen(model.numElements,model.R,model.hubRadius,model.preSweep,model.preCurve,model.edgeEAOffset,model.flapEAOffset,model.elementOrder,model.mesh);
model.numNodes = length(model.mesh.x);
model.numDofPerNode = 6;
model.totalNumDof = model.numNodes * model.numDofPerNode;
%------------------------------
t_mesh = toc;
disp('Elapsed time for mesh generation(s):');
disp(t_mesh);

%-------- ANALYSIS ------------
   displacement = zeros(model.totalNumDof,1);
   OmegaStart = 0.0;
   %Do nonlinear iteration if needed
    if(strcmp(analysisType,'F'))    
        [displacement,staticAnalysisSuccessful]=staticAnalysis(model,displacement,Omega,OmegaStart,smOmega,0,analysisType);
    elseif(strcmp(analysisType,'P') | strcmp(analysisType,'D')) 
        staticAnalysisSuccessful = true; 
    end
   if(~staticAnalysisSuccessful)
       error('Static non-linear analysis not successful. Maximum iterations exceeded.');
   end
   
   %perform modal analysis
   [freq,damp,phase1,phase2,imagCompSign,convergedCheck]=linearAnalysis(model,displacement,Omega,smOmega,aeroLoadsFlag,analysisType); 
%------------------------------

   %write output
   writeOutput(freq,damp,phase1,phase2,imagCompSign,model.R,convergedCheck,fid,analysisType,1);

fclose(fid);
end