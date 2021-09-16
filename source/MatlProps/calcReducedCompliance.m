function y = calcReducedCompliance(E1,E2,NU12,G12)
%ReducedCompliance This function returns the reduced compliance
% matrix for fiber-reinforced materials.
% There are four arguments representing four
% material constants. The size of the reduced
% compliance matrix is 3 x 3.
%
% See Chapter 4 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information
y = [1/E1    -NU12/E1 0 ; 
    -NU12/E1  1/E2    0 ;
     0        0       1/G12];
end