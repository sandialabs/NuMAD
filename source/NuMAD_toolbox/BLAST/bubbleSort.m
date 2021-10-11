function [A,origMap,posIndex] = bubbleSort(A)
%bubbleSort Sorts the vector A in ascending order
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [A,origMap,posIndex] = bubbleSort(A)
%                    
%   This function accepts a vector A, sorts the vector in ascending order,
%   outputting the sorted vector, a map to the original ordering, and the
%   index at which a positive value first occurs.
%
%      input:
%      A        = vector to be sorted
%
%      output:
%      A        = sorted vector
%      origMap  = map of sorted vector to original ordering
%      posIndex = index at which positive value first occurs

Aorig=A;
    len = length(A);
    swapped = true;
    origMap = 1:len;
    while(swapped)
    swapped = false;

    for i=1:len-1
        if(A(i+1) < A(i))
            temp = A(i);
            A(i) = A(i+1);
            A(i+1) = temp;
            swapped = true;
            
            temp2 = origMap(i);
            origMap(i) = origMap(i+1);
            origMap(i+1) = temp2;
        end
    end
    end
    
%     posIndex = length(A)/2+1;
%     [posIndex] = findPositiveCrossOver(A);
posIndex = 1;
end


function [posIndex] = findPositiveCrossOver(A)
%This function finds the first insance of a positive entry in a vector
%      input:
%      A        = input vector
%
%      output:
%      posIndex = index of first instance of a positive entry in A
    for i=1:length(A)
        if(A(i) >= 0)
            posIndex = i;
            break;
        end
        
    end
    
end
