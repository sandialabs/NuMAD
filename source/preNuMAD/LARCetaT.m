function[etaT]=LARCetaT(alp0)
%"In the absence of biaxial test data, ?L can be estimated from the longitudinal
%and transverse shear strengths." Failure Criteria for FRP Laminates in Plane Stress
% Carlos G. Dávila, Pedro P. Camanho

% Compute the coefficient of transverse influence required for Larc failure criteria.

%INPUTS
% alp0 - Material fracture angle, degrees
etaT=-1/tand(2*alp0);