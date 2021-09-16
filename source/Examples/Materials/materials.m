% Example file for calculation of effective laminate properties for a
% stack of individual laminae.
%
% Inputs: set layer orientations and total laminate thickness (this example
% assumes each layer is of equal thickness).  Fiber and matrix densities
% and fiber volume fraction. Calculate E1, E2, nu12, G12 of laminae based 
% either on micromechanics or given values.  The script does the rest!
%
% Output to workspace: Ex, Ey, NUxy, and Gxy for the entire
% laminate based on A matrix (assumption is "balanced symmetric" layup)

clear all
theta=[-45 45 zeros(1,23) 45 -45];  % theta for each layer
H=2.82;  % Total laminate thickness, mm
rhof=1800; % kg/m^3 or 1000*g/cc or 1000*specific_gravity
rhom=1200; % kg/m^3 or 1000*g/cc or 1000*specific_gravity
Vf=0.55;

if 0  % Micromechanics example
    E1f  = 233e9;
    E2f  = 14.8e9;
    G12f = 8.96e9;
    G23f = 8.27e9;
    Em   = 3.45e9;
    num  = 0.36;
    nu12f= 0.200;
    nu23f= 0.200;
    Gm=Em/(2*(1+num));
    
    disp('Calculated laminae properties:')
    E1=calcE1(Vf,E1f,Em);
    E2=calcE2(Vf,E2f,Em,0,0,0,num,0,1,5);
    nu12=calcNU12(Vf,nu12f,num);
    G12=calcG12(Vf,G12f,Gm,0,1);
else  % Known laminae properties example
    % Ref[1]: Isaac M. Daniel and Ori Ishai.  Engineering Mechanics of Composite Materials. Oxford University Press, 1994.
    disp('Given laminae properties:')
    E1 =  114.5e9;  % tweaked to match Ex of 100.1 GPa
    E2 = E1/13.64;  % from Ref[1], Table 2.3
    nu12 = 0.27;  % from Ref[1], Table 2.6
    G12 = E1/19.1;  % from Ref[1], Table 2.3
end
disp(sprintf('E1=%6.2fGPa, E2=%6.2fGPa, nu12=%4.3f, G12=%6.2fGPa\n',E1/1e9,E2/1e9,nu12,G12/1e9))

% calculate A,B,D for laminate
N=length(theta)+1;
z=linspace(-H/2,H/2,N);  % layer boundaries, distances from midplane

S=calcReducedCompliance(E1,E2,nu12,G12);  % "reduced" = assumption of plane stress
Q=inv(S); % calculate reduced stiffness matrix
A=zeros(3,3);  % initialize A
B=zeros(3,3);  % initialize B
D=zeros(3,3);  % initialize D
for i=1:length(theta) % Loop through each layer in the laminate
    T=calcT(theta(i));
    Qbar=inv(T)*Q*T;
    A=calcAmatrix(A,Qbar,z(i),z(i+1));
    B=calcBmatrix(B,Qbar,z(i),z(i+1));
    D=calcDmatrix(D,Qbar,z(i),z(i+1));
end

% Calculate effective laminate properties from A matrix
% Assumes "balanced symmetric" layup
disp('=>Calculated laminate properties:')
rho=Vf*rhof+(1-Vf)*rhom;
Ebarx=calcEbarx(A,H);
Ebary=calcEbary(A,H);
nubarxy=calcNUbarxy(A,H);
nubaryx=calcNUbaryx(A,H);
Gbarxy=calcGbarxy(A,H);
disp(sprintf('rho=%4.0fkg/m^3, E1=%6.2fGPa, E2=%6.2fGPa, nu12=%4.3f, G12=%6.2fGPa\n',rho,Ebarx/1e9,Ebary/1e9,nubarxy,Gbarxy/1e9))

