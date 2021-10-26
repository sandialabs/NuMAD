function [freqSorted,dampSorted,phase1Sorted,phase2Sorted,imagCompSignSorted,convergedCheckSorted] = sortModes(freq,damp,phase1,phase2,imagCompSign,convergedCheck)
%sortModes sorts modes in ascending frequency
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [freqSorted,dampSorted,phase1Sorted,phase2Sorted,imagCompSignSorted,
%    convergedCheckSorted] = sortModes(freq,damp,phase1,phase2,
%                                      imagComponentSign,convergedCheck)
%
%   This function sorts modes by ascending freqency. Frequency, damping
%   ratio, 0 deg phase mode shape, 90 deg phase mode shapes, imaginary sign
%   component, and converged frequency flag array are sorted with respect
%   to ascending frequency.
%
%      input:
%      freq           = array of unsorted frqeuencies
%      damp           = array of unsorted damping ratios
%      phase1         = array of unsorted 0 degree phase mode shapes
%      phase2         = array of unsorted 90 degree phase mode shapes
%      imagCompSign   = array of unsorted imaginary component 
%      convergedCheck = array of unsorted converged frequency check flags
%      
%      output:
%      freqSorted            = array of sorted frequencies
%      dampSorted            = array of sorted damping ratios
%      phase1Sorted          = array of sorted 0 degree phase mode shapes
%      phase2Sorted          = array of sorted 90 degree phase mode shapes
%      imagCompSignSorted    = array of sorted imaginary component signs
%      convergedCheckSorted  = array of converged frequency check flags

[freq,map,posIndex] = bubbleSort(freq);

    for i=posIndex:1:posIndex+19
        dampSorted(i) = damp(map(i));
        freqSorted(i) = freq(i);
        phase1Sorted(:,:,i) = phase1(:,:,map(i));
        phase2Sorted(:,:,i) = phase2(:,:,map(i));
        imagCompSignSorted(i) = imagCompSign(map(i));
        convergedCheckSorted(i) = convergedCheck(map(i));
    end


end
