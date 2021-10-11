function y = calcDmatrix(D,Qbar,z1,z2)
%Dmatrix This function returns the [D] matrix
% after the layer k with stiffness [Qbar]
% is assembled.
% D - [D] matrix after layer k
% is assembled.
% Qbar - [Qbar] matrix for layer k
% z1 - z(k-1) for layer k
% z2 - z(k) for layer k
for i = 1 : 3
    for j = 1 : 3
        D(i,j) = D(i,j) + Qbar(i,j)*(z2^3 -z1^3)/3;
    end
end
y = D;
end