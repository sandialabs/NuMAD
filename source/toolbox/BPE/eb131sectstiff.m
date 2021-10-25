function [storeData] = eb131sectstiff(nodeData,storeData)
% sandia equiv beam study
% calculates the 6x6 stiffnesses of each section (at center of element)
% 27  nov 2001, corrected L(4,2): 6=> 3.
% latest: 10 Sept 2002
% 20 nov 2003, v27 for use with Legacy Numad
%  14 april 2004, outer profile added
%  22 july 2005, reference to profile removed (solid elements)
%  02 september 2009, brr, cleaned up for use in NuMAD
%  08 dec 2010, changed variable names to match what's in 2006 W.E. paper. - b.resor
%

nbeamnode = nodeData.nbeamnode;
zbeamnode = nodeData.zbeamnode;

storestiffsym = storeData.storestiffsym;

storesectstiff = zeros((nbeamnode-1)*8,6);
% create common E matrix
E = zeros(6,6);
E(1,5) = 1;
E(2,4) = -1;

for ielement = 1:nbeamnode-1       % loop on beam elements
    
    length = zbeamnode(ielement+1)-zbeamnode(ielement);
    vector = [1 1 1 1 1 1]*length;
    % create H matrix
    H = diag(vector);
    H(4,2) = -(length^2)/2;
    H(5,1) = -H(4,2);
    % create L matrix
    Q = diag(vector)*length/2;
    Q(4,2) = -(length^3)/3;
    Q(5,1) = -Q(4,2);
    % inputs for Lyapunov eqtn: -C = AX + XB
    B = H*inv(Q);
    A = E;
    i1 = (ielement-1)*8;
    C = storestiffsym(i1+1:i1+6,1:6);
    % C = storestiff(i1+1:i1+6,1:6);
    C = inv(C)*inv(Q);
    X = lyap(A,B,C);  %  solve Lyapunov equation for X =-inv(stiffness matrix)
    matrix = -inv(X);
    storesectstiff(i1+1:i1+6,1:6) = matrix;
    
    %  brr added the following line
    storesectflex(i1+1:i1+6,1:6) = X;
    
end  % end of element loop

% save section stiffnesses for later use
save sectstiffmatrix.txt storesectstiff -ascii
save sectflexmatrix.txt storesectflex -ascii

storeData.storesectstiff = storesectstiff;

end
