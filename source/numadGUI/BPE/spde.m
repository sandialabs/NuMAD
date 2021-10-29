function X = spde(A,B,showCond)
% This function solves an inverse problem assuming the that resulting
%   matrix should be symmetric positive definite.  Minimizes the area of
%   the error.  The equation is solved using the symmetric positive 
%   definite esimation problem method. Look in the reference for more 
%   information.
%
% INPUTS:
%   A - The data matrix in AX=B.
%   B - The data matrix in AX=B.
%   showCond - A boolean that if true will show the condition number of the
%       data matrix A.
%
% OUTPUTS:
%   X - the solution matrix in AX=B.
%
% This is based off Yixin Chen's thesis
% http://cs.olemiss.edu/~ychen/publications/others/thesis.pdf
  
if nargin > 2 && showCond
    dispCond = cond(A)
end
P = A'*A;
Q = B'*B;
[U,D1] = eig(P);
s1 = sqrt(D1);
tmp = U*s1;
[V,D2] = eig(tmp'*Q*tmp);
s2 = sqrt(D2);

% Take the inverse of s1
iS1 = diag(1./diag(s1));
tmp = U*iS1*V;
X = tmp*s2*tmp';