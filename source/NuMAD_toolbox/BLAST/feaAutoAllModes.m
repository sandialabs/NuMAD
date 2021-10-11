function [convergedFreq,convergedDamp]=feaAutoAllModes(inStruct,tol,analysisType)
%feaAutoAllModes Performs automated flutter analysis
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [convergedFreq,convergedDamp]=feaAutoAllModes(infile,outfile,
%                                 OmegaArray,tol,analysisType)
%
%   This function reads input, generates a beam element mesh, performs
%   static "spin up" analysis, and modal analysis of system. Flutter
%   analysis is conducted in an automated manner for the first few modes of
%   a system. A "p-k" iteration process is utilized until convergence is
%   reached for a particular mode.
%
%      input:
%      inStruct      = struct containing input filename
%      tol           = tolerance to assess convergence of a mode
%      analysisType  = char denoting stability analysis type (F-rotating
%                      flutter analysis, P-parked flutter analysis, 
%                      D-rotating divergence analysis)
%      
%      output:
%      convergedFreq = array of converged frequencies
%      convergedDamp = array of converged damping ratios

%set analysis type
if(strcmp(analysisType,'F'))
    disp('Running in rotating flutter mode');
elseif(strcmp(analysisType,'P'))
    disp('Running in parked flutter mode');
elseif(strcmp(analysisType,'D'))
    disp('Redirecting to fea() for divergence analysis ...');
    fea(inStruct,0,0,1,analysisType);
    return;
else
    error('Analysis type not recognized. Exiting.');
end

%set analysis flags and input parameters
aeroLoadsFlag = 2;  %1=Lobitz complex aero matrices,  2=Wright and Cooper real valued representation
smOmega = 1.0;
OmegaArray = inStruct.OmegaArray;
tic

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

%Check FAST and NuMAD data for blade length  consistency
if(abs(model.bladeLength-inStruct.numadBladeLen)>1.0e-3)
    disp('WARNING: Inconsistency in NuMAD and FAST Blade Lengths.');
end

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


displacement = zeros(model.totalNumDof,1);
OmegaStart = 0.0;

%loop over operation condition array (windspeed, rotor speed, etc)
for i=1:length(OmegaArray)

%this step sets up an initial solution to iterate from for converged
%frequencies

%-------- ANALYSIS ------------
%Do nonlinear iteration if needed
if(strcmp(analysisType,'F'))
    disp(['Running at rotor speed ',num2str(OmegaArray(i)),' RPM']);
   [displacement,staticAnalysisSuccessful]=staticAnalysis(model,displacement,OmegaArray(i),OmegaStart,smOmega,0,analysisType);
elseif(strcmp(analysisType,'P')) 
   staticAnalysisSuccessful = true;
   disp(['Running at freestream velocity ',num2str(OmegaArray(i))]);
end

   if(staticAnalysisSuccessful)
   %perform modal analysis
   [freqOrig,damp,phase1,phase2,imagCompSign,convergedCheck]=linearAnalysis(model,displacement,OmegaArray(i),smOmega,aeroLoadsFlag,analysisType); 
    %write output
%    writeOutput(freqOrig,damp,phase1,phase2,imagCompSign,model.R,convergedCheck,fid,analysisType,0);
   OmegaStart = OmegaArray(i);

   %loops over modes of interest (first 10 modes) to converged frequency at
   %given operation condition
for j=1:length(freqOrig)   
   converged = 0;
   if(OmegaArray(i) == 0.0)
       converged = 1;
       convergedFreq(i,:) = freqOrig;
       convergedDamp(i,:) = damp;
   end
   smOmega = freqOrig(j);
   while(~converged)
   %perform modal analysis
   [freq,damp,phase1,phase2,imagCompSign,convergedCheck]=linearAnalysis(model,displacement,OmegaArray(i),smOmega,aeroLoadsFlag,analysisType); 
    %write output
%    writeOutput(freq,damp,phase1,phase2,imagCompSign,model.R,convergedCheck,fid,analysisType,0);
%------------------------------
    
   %new guess for smOmega
   if(smOmega-freq(j)<tol)
       converged = true;
       convergedFreq(i,j) = freq(j);
       convergedDamp(i,j) = damp(j);
       convergedImagSign(i,j) = imagCompSign(j);
       if(damp(j)<0)
           fprintf('    * Negative damping present\n');
       end
       break;
   else
   smOmega = 0.5*(freq(j)+smOmega);
   end
   
   end
end

   end
   
   %option to output available data if nonlinear static analysis fails to
   %converge
   if(~staticAnalysisSuccessful)
       lastSuccessfulIndex = i-1;
       break;
   else
       lastSuccessfulIndex = i;
   end
end

%plot frequency and damping vs. omega/freestream velocity
omegaPlot = OmegaArray(1:lastSuccessfulIndex);
plotFreqDampCurves(omegaPlot,convergedDamp,convergedFreq,analysisType);


%print output summary file header
outfile = removeFileExtensionFromFilename(inStruct.outFile);
outfid = printOutputSummaryHeaderStruct(outfile,inStruct,analysisType);

%search for negative damping and interpolate for "cross-over" operating
%condition
identifyPotentialInstabilities(convergedDamp,convergedFreq,omegaPlot,...
                               convergedImagSign,analysisType,outfid);
                           
%save frequency and damping information to .mat file
save(outfile,'omegaPlot','convergedDamp','convergedFreq','analysisType');


end

function [newoutfile] = removeFileExtensionFromFilename(outfile)

    len = length(outfile);
    
    index = len+1;
    for i=1:len
       if(strcmp(outfile(i),'.'))
          index = i;
          break
       end
    end
    
    newoutfile = outfile(1:index-1);
  
end