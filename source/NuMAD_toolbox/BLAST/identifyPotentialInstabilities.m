function identifyPotentialInstabilities(damp,freq,omega,imagCompSign,analysisType,fid)
%identifyPotentialInstabilities  identifies potential instablilities
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   identifyPotentialInstabilities(damp,freq,omega,analysisType,fid)
%                    
%   This funciton searches for potential instabilities by search for
%   negative damping conditions. A prediction of the instability onset is
%   calculated by interpolating to find the operating condition associated
%   with the zero-damping value at the cross-over point.
%
%      input:
%      damp             = array of damping ratio of various modes over 
%                         operating conditions
%      freq             = array of frequency of various modes over
%                         operating conditions
%      omeg             = array of operating conditions
%      analysisType     = character denoting analysis type
%      fid              = file id for output summary file (.sum)
%
%      output:
%      none

[numRPM,numModes] = size(freq);
    index = 1;
    
    instabRpm = [];
    
for i=1:numModes
    negIndex = 0;
   for j=1:numRPM
      if(damp(j,i) < 0 && omega(j) > 0 && abs(damp(j,i))>1.0e-8)
         negIndex = j;
         break;
      end
   end
   
   if(negIndex >1)
       damp1 = damp(negIndex-1,i);
       damp2 = damp(negIndex,i);
       
       if(abs(damp1)<1.0e-8)
            damp1= 0;
       end
       
       freq1 = freq(negIndex-1,i);
       freq2 = freq(negIndex,i);
       
       rpm1 =  omega(negIndex-1);
       rpm2 =  omega(negIndex);
       
       m = (damp2-damp1)/(rpm2-rpm1);
       
       rpm = rpm1 - damp1/m;
       
       freqI = interp1([rpm1 rpm2],[freq1 freq2],rpm);
       
       instabRpm(index) = rpm;
       instabFreq(index) = freqI;
       instabMode(index) = i;
       instabSign(index) = imagCompSign(negIndex,i);
       index = index + 1;
   elseif(negIndex == 1 && omega(1)~=0.0)
       instabRpm(index) = omega(1);
       instabFreq(index) = freq(1,i);
       instabMode(index) = i;
       instabSign(index) = imagCompSign(negIndex,i);
       index = index + 1;
   end
   
end

[rpms,map] = sort(instabRpm);

fprintf('\n\n**************************************************\nANALYSIS SUMMARY:\n\n');
fprintf(fid,'\n\n*************************************************\nANALYSIS SUMMARY:\n\n');
for i=1:length(rpms)
    
   fprintf('Potential instability identified at %4.2f',rpms(i)); 
   fprintf(fid,'Potential instability identified at %4.2f',rpms(i));
   if(strcmp(analysisType,'F'))
      fprintf(' RPM.\n');
      fprintf(fid,' RPM.\n');
   end
   if(strcmp(analysisType,'P'))
      fprintf(' m/s.\n');
      fprintf(fid,' m/s.\n');
   end
   fprintf('Mode %i:\t Frequency %4.2f (Hz)\n',instabMode(map(i)),instabFreq(map(i)));
   fprintf(fid,'Mode %i:\t Frequency %4.2f (Hz)\n',instabMode(map(i)),instabFreq(map(i)));
   if(instabSign(map(i))<0)
        fprintf('\t\t~~Negative Frequency~~');
        fprintf(fid,'\t\t~~Negative Frequency~~');
   end
       fprintf('\n\n');
       fprintf(fid,'\n\n');

end

if(isempty(instabRpm))
    fprintf('No potential instabilities found.\n');
    fprintf(fid,'No potential instabilities found.\n');
end

fclose(fid);

end

