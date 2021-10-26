% References
% [1] 88-12-15/AB Rotor Blade.pdf
% [2] Ever J. Barbero. Introduction to Composite Design. Taylor and
%     Francis, 1999.
% [3] 88-12-15/AB Rotor Blade.pdf (Encl.5.A4)
%

clear all

t1200=0.078*12;  % mm, thickness of UD1200

%% define Micromechanics properties; fiber and matrix
Vf   = 0.50;  % Page 1 [3]
rhof = 2550; % kg/m^3 or 1000*g/cc or 1000*specific_gravity; page 5 [1]
rhom = 1200; % kg/m^3 or 1000*g/cc or 1000*specific_gravity; a guesstimate; back calcuated from [3]
% Oriented fiber material properties
E1f  = 70e9;  % from p.5 [1]
E2f  = E1f;  % Isotropic assumption
nuf  = 0.22; % from Table 2.1 [2]
nu12f= nuf;
nu23f= nuf;  % assumption
Gf   = E1f/(2*(1+nuf)); % Isotropic assumption
G12f = Gf;
G23f = Gf; % assumption
Em   = 4e9;  % back calculated from [3]
num  = 0.38;  % polyester resin, Table 2.4 [2]
Gm   = Em/(2*(1+num));  % Isotropic assumption
Eta  = 0;
nu21f= 0;

%% M300 Random Mat
t_m300=300/1200*t1200;
Ex_m300=8.8e9;
Ey_m300=Ex_m300;
nuxy_m300=0.3;
Gxy_m300=Ex_m300/(2*(1+nuxy_m300));
rho_m300=1200*.5+2550*.5;
disp('M300 -- given laminae properties (from Vestas):')
m300=[t_m300,rho_m300,Ex_m300/1e9,Ey_m300/1e9,nuxy_m300,Gxy_m300/1e9];
disp(sprintf('t=%5.3fmm, rho=%4.0fkg/m^3, Ex=%6.2fGPa, Ey=%6.2fGPa, nuxy=%4.3f, Gxy=%6.2fGPa\n',m300))

%% M100/UD700
%M100
t_m100=100/1200*t1200;
Ex_m100=Ex_m300;
Ey_m100=Ex_m100;
nuxy_m100=nuxy_m300;
Gxy_m100=Gxy_m300;
rho_m100=rho_m300;
disp('M100 -- given laminae properties (based on Vestas):')
m100=[t_m100,rho_m100,Ex_m100/1e9,Ey_m100/1e9,nuxy_m100,Gxy_m100/1e9];
disp(sprintf('t=%5.3fmm, rho=%4.0fkg/m^3, E1=%6.2fGPa, E2=%6.2fGPa, nu12=%4.3f, G12=%6.2fGPa\n',m100))

%UD700
t_ud700=700/1200*t1200;
E1_ud700=calcE1(Vf,E1f,Em);
E2_ud700=calcE2(Vf,E2f,Em,Eta,nu12f,nu21f,num,E1f,1,4);
nu12_ud700=calcNU12(Vf,nu12f,num);
G12_ud700=calcG12(Vf,G12f,Gm,0,3);
rho_ud700=Vf*rhof+(1-Vf)*rhom;
disp('UD700 -- calculated laminae properties (from CLT):')
ud700=[t_ud700,rho_ud700,E1_ud700/1e9,E2_ud700/1e9,nu12_ud700,G12_ud700/1e9];
disp(sprintf('t=%5.3fmm, rho=%4.0fkg/m^3, E1=%6.2fGPa, E2=%6.2fGPa, nu12=%4.3f, G12=%6.2fGPa\n',ud700))

theta=[0 0];  % theta for each layer
H=t_m100+t_ud700;  % Total laminate thickness, mm

% calculate A,B,D for laminate
z=[-H/2 H/2-t_ud700 H/2];  % layer boundaries, distances from midplane
A=zeros(3,3);  % initialize A
B=zeros(3,3);  % initialize B
D=zeros(3,3);  % initialize D

i=1;  % Layer 1
S=calcReducedCompliance(Ex_m100,Ey_m100,nuxy_m100,Gxy_m100);  % "reduced" = assumption of plane stress
Q=inv(S); % calculate reduced stiffness matrix
T=calcT(theta(i));
Qbar=inv(T)*Q*T;
A=calcAmatrix(A,Qbar,z(i),z(i+1));
B=calcBmatrix(B,Qbar,z(i),z(i+1));
D=calcDmatrix(D,Qbar,z(i),z(i+1));

i=2;  % Layer 2
S=calcReducedCompliance(E1_ud700,E2_ud700,nu12_ud700,G12_ud700);  % "reduced" = assumption of plane stress
Q=inv(S); % calculate reduced stiffness matrix
T=calcT(theta(i));
Qbar=inv(T)*Q*T;
A=calcAmatrix(A,Qbar,z(i),z(i+1));
B=calcBmatrix(B,Qbar,z(i),z(i+1));
D=calcDmatrix(D,Qbar,z(i),z(i+1));

% Calculate effective laminate properties from A matrix
% Assumes "balanced symmetric" layup
disp('=>M100/UD700 -- Calculated laminate properties (based on CLT):')
rho=t_m100/H*rho_m100+t_ud700/H*rho_ud700;
Ebarx=calcEbarx(A,H);
Ebary=calcEbary(A,H);
nubarxy=calcNUbarxy(A,H);
nubaryx=calcNUbaryx(A,H);
Gbarxy=calcGbarxy(A,H);
disp(sprintf('t=%5.3fmm, rho=%4.0fkg/m^3, Ex=%6.2fGPa, Ey=%6.2fGPa, nuxy=%4.3f, Gxy=%6.2fGPa\n',...
    H,rho,Ebarx/1e9,Ebary/1e9,nubarxy,Gbarxy/1e9))

%% DB600/UD300 Triax
t_ud300=300/1200*t1200;
E1_ud300=E1_ud700;
E2_ud300=E2_ud700;
nu12_ud300=nu12_ud700;
G12_ud300=G12_ud700;

theta=[-45 0 45];  % theta for each layer
H=t_ud300*3;  % Total laminate thickness, mm

% calculate A,B,D for laminate
z=[-H/2  H/2-t_ud300*2 H/2-t_ud300 H/2];  % layer boundaries, distances from midplane
A=zeros(3,3);  % initialize A
B=zeros(3,3);  % initialize B
D=zeros(3,3);  % initialize D

i=1;  % Layer 1
S=calcReducedCompliance(E1_ud300,E2_ud300,nu12_ud300,G12_ud300);  % "reduced" = assumption of plane stress
Q=inv(S); % calculate reduced stiffness matrix
T=calcT(theta(i));
Qbar=inv(T)*Q*T;
A=calcAmatrix(A,Qbar,z(i),z(i+1));
B=calcBmatrix(B,Qbar,z(i),z(i+1));
D=calcDmatrix(D,Qbar,z(i),z(i+1));

i=2;  % Layer 2
S=calcReducedCompliance(E1_ud300,E2_ud300,nu12_ud300,G12_ud300);  % "reduced" = assumption of plane stress
Q=inv(S); % calculate reduced stiffness matrix
T=calcT(theta(i));
Qbar=inv(T)*Q*T;
A=calcAmatrix(A,Qbar,z(i),z(i+1));
B=calcBmatrix(B,Qbar,z(i),z(i+1));
D=calcDmatrix(D,Qbar,z(i),z(i+1));

% Calculate effective laminate properties from A matrix
% Assumes "balanced symmetric" layup
disp('=>DB600/UD300 Triax -- Calculated laminate properties (based on CLT):')
rho=rho_ud700;
Ebarx=calcEbarx(A,H);
Ebary=calcEbary(A,H);
nubarxy=calcNUbarxy(A,H);
nubaryx=calcNUbaryx(A,H);
Gbarxy=calcGbarxy(A,H);
disp(sprintf('t=%5.3fmm, rho=%4.0fkg/m^3, Ex=%6.2fGPa, Ey=%6.2fGPa, nuxy=%4.3f, Gxy=%6.2fGPa\n',...
    H,rho,Ebarx/1e9,Ebary/1e9,nubarxy,Gbarxy/1e9))
