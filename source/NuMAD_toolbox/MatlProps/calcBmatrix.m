function y = calcBmatrix(B,Qbar,z1,z2)
%Bmatrix This function returns the [B] matrix
% after the layer k with stiffness [Qbar]
% is assembled.
% B - [B] matrix after layer k
% is assembled.
% Qbar - [Qbar] matrix for layer k
% z1 - z(k-1) for layer k
% z2 - z(k) for layer k
for i = 1 : 3
    for j = 1 : 3
        B(i,j) = B(i,j) + Qbar(i,j)*(z2^2 -z1^2)/2;
    end
end
y = B;
end