function X = spp(A,B,showCond)
% This function solves an inverse problem assuming the that resulting
%   matrix should be symmetric positive definite.  This solution technique
%   is based off the symmetric procrustes problem which solves a modified
%   version of the Lyapunov equation.
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
X = lyap(A'*A,-1*(A'*B+B'*A));