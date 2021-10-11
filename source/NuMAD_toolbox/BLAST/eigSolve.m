function [eigVec,eigVal] = eigSolve(M,C,K,flag)
%eigSolve Performs eigensolve on discretized equations of motion
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [eigVec,eigVal] = eigSolve(M,C,K,flag)
%                    
%   This function accepts mass, damping, and stiffness matrices of the
%   system, calculates the state-space representation, and performs an
%   eigensolve, outputting eigenvectors and eigenvalues.
%
%      input:
%      M        = system mass matrix
%      C        = system damping matrix
%      K        = system stiffness matrix
%
%      output:
%      eigVec   = eigenvectors of the state-space system representation
%      eigVal   = eigenvalues of the state-space system representation


	[len,dum] = size(M);
	sysMat = [zeros(len), eye(len);
			  -(M^-1)*K, -(M^-1)*C];
          
     	
    if(flag==1)
        [eigVec,eigVal] = eig(sysMat);		  
    end
    if(flag==2)
         sysMat=sparse(sysMat);
        [eigVec,eigVal] = eigs(sysMat,20,'SM');		  
    end
    
end