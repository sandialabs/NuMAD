function y = calcOrthotropicCompliance(E1,E2,E3,NU12,NU23,NU13,G12,G23,G13)
%OrthotropicCompliance This function returns the compliance matrix
% for orthotropic materials. There are nine
% arguments representing the nine independent
% material constants. The size of the compliance
% matrix is 6 x 6.
%
% See Chapter 2 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

y = [ 1/E1   -NU12/E1  -NU13/E1 0     0     0 ;
    -NU12/E1  1/E2     -NU23/E2 0     0     0 ;
    -NU13/E1 -NU23/E2   1/E3    0     0     0 ;
     0        0         0       1/G23 0     0 ;
     0        0         0       0     1/G13 0 ;
     0        0         0       0     0     1/G12];
end