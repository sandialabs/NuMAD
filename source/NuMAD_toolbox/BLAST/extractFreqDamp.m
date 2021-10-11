function [freq,damp,phase1,phase2] = extractFreqDamp(val,vec,numDOFPerNode,dofVector,analysisType)
%ExtractFreqDamp Extracts frequency and damping info from eigenvalues
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [freq,damp,phase1,phase2] = extractFreqDamp(val,vec,numDOFPerNode,
%                               dofVector,analysisType)
%                    
%   This function extracts frequency, damping, and phase information from
%   the eigenvalues. 
%
%      input:
%      val           = eigenvalue
%      vec           = eigenvector
%      numDOFPerNode = number of degrees of freedom per node
%      dofVector     = degree of freedom vector associated with vec
%      analysisType  = char for definiing stability analysis type
%
%      output:
%      freq           = frequency of mode (Hz)
%      damp           = damping ratio for mode
%      phase1         = 0 degree (in-phase) mode shape
%      phase2         = 90 degree (out-of-phase) mode shape
     
     freq = abs(imag(val))/(2*pi);
     damp = -2.0*real(val)/abs(imag(val));
     
     if(abs(imag(val)) < 1.0e-6)
         freq = sqrt(abs(real(val)))/(2*pi);
         damp = 0.0;
     end
     
     if(strcmpi(analysisType,'D'))
         if(isreal(val))
             freq = val;
         else
             freq = abs(imag(val));
         end
         freq = freq/(2*pi)*60;
     end
     
     
  
     [len,numModeShapes] = size(vec);
      sortedModes=zeros(len/(numDOFPerNode*2)+1,numDOFPerNode);

    for i=1:len/2
       sortedModes(dofVector(i,2),dofVector(i,3)) = vec(i);
    end
     
     phase1 = real(sortedModes);
     phase2 = imag(sortedModes);
     
     max1=max(max(abs(phase1)));
     max2=max(max(abs(phase2)));
     
     if(max1>max2)
         maxOverall = max1;
     else
         maxOverall = max2;
     end
     
     phase1 = phase1./maxOverall;
     phase2 = phase2./maxOverall;
     
     if(abs(min(min(phase1))+1)<1.0e-4 || abs(min(min(phase2))+1)<1.0e-4)
         phase1 = -1*phase1;
         phase2 = -1*phase2;
         
     end
     

end