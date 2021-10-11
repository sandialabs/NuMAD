function[etaL]=LARCetaL(SL,YC,alp0)
%"In the absence of biaxial test data, ?L can be estimated from the longitudinal
%and transverse shear strengths." Failure Criteria for FRP Laminates in Plane Stress
% Carlos G. Dávila, Pedro P. Camanho

% Compute the coefficient of longitudinal influence required for Larc failure criteria.

%INPUTS
%   SL - Lateral shear strength 
%   YC - Transverse compressive strength 
% alp0 - Material fracture angle, degrees

etaL=-SL*cosd(2*alp0)/(YC*cosd(alp0)^2);