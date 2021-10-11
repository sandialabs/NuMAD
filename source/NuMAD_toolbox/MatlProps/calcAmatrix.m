function y = calcAmatrix(A,Qbar,z1,z2)
%Amatrix This function returns the [A] matrix
% after the layer k with stiffness [Qbar]
% is assembled.
% A - [A] matrix after layer k
% is assembled.
% Qbar - [Qbar] matrix for layer k
% z1 - z(k-1) for layer k
% z2 - z(k) for layer k
for i = 1 : 3
    for j = 1 : 3
        A(i,j) = A(i,j) + Qbar(i,j)*(z2-z1);
    end
end
y = A;
end