function [lambda] = calculateLambda(theta1,theta2,theta3 )
%calculateLambda Calculates transformation matrix from element to hub frame
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [lambda] = calculateLambda(theta1,theta2,theta3 )
%                    
%   This function calculates a transformation matrix to transform the
%   element degree of freedom vector (12 DOFs) from the element frame to
%   the hub frame. The transformation matrix is constructed via the
%   direction cosine matrices of a 3-2-1 Euler rotation sequence.
%
%      input:
%      theta1        = angle (rad) of rotation for "1" of 3-2-1 sequence
%      theta2        = angle (rad) of rotation for "2" of 3-2-1 sequence
%      theta3        = angle (rad) of rotation for "3" of 3-2-1 sequence

%      output:
%      lambda        = 12 x 12 transformation matrix

    dcm1 = [1 0 0;
            0 cos(theta1) sin(theta1);
            0 -sin(theta1) cos(theta1)];
    
    dcm2 = [cos(theta2) 0 -sin(theta2);
            0 1 0;
            sin(theta2) 0 cos(theta2)];
    
    dcm3 = [cos(theta3) sin(theta3) 0
            -sin(theta3) cos(theta3) 0
             0           0           1];
     
         % 3 (sweep) - 2(dihedral) - 1(twist) Euler Sequence
     dcm = dcm1*dcm2*dcm3;
     
     zm = zeros(3);
     
     lambda = [dcm,zm,zm,zm;
                zm,dcm,zm,zm;
                zm,zm,dcm,zm;
                zm,zm,zm,dcm];

end

