function [phase1,phase2] = readResults(resultsFile,numModes,numNodes)
%readResults Reads results file and outputs in/out of phase mode shapes
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [phase1,phase2] = readResults(resultsFile,numModes,numNodes)
%                    
%   This function reads the BLAST results and outputs mode shape
%   information for in-phase and out-of-phase shapes.
%
%      input:
%      resultsFile    = BLAST results file
%      numModes       = number of modes in results file
%      numNodes       = number of nodes in mesh
 
%      output:
%      phase1             = in-phase mode shapes
%      phase2             = out-of-phase mode shapes


fid = fopen(resultsFile,'r');

modecount = 1;

% while(~feof(fid))
while(modecount <=numModes)
dum = fscanf(fid,'%s',17);
for j=1:numNodes
   temp1 = fscanf(fid,'%e',6);
   phase1(j,:,modecount) = temp1';
end

dum = fscanf(fid,'%s',10);
for j=1:numNodes
   temp2 = fscanf(fid,'%e',6);
   phase2(j,:,modecount) = temp2';
end
modecount = modecount + 1;
end
fclose(fid);

end

