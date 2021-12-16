function [Ke,Me,Ce,Fe] = calculateEulerBeamElement(elementOrder,x,y,z,xloc,hubRadius,disp,sectionProps,sweepAngle,coneAngle,aeroSweepAngle,Omega,smOmega,airDensity,aeroLoadsFlag,modalFlag,analysisType)
%calculateEulerBeamElement Calculates element matrices and vector
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [Ke,Me,Ce,Fe] = calculateEulerBeamElement(elementOrder,x,y,z,xloc,disp,
%                   sectionProps,sweepAngle,aeroSweepAngle,Omega,smOmega,
%                   airDensity,aeroLoadsFlag,modalFlag,analysisType)
%
%   This function calculates the element stiffness, mass, and damping
%   matrices and load vector for an Euler-Bernoulli beam element with
%   rotational effects.
%
%      input:
%      elementOrder     = order of element (1-linear currently implemented)
%      x                = x coordinates of element nodes
%      y                = y coordinates of element nodes
%      z                = z coordinates of element nodes
%      xloc             = local element coordinate of element nodes
%      hubRadius        = radius of hub
%      disp             = displacement of element nodes
%      sectionProps     = section properties for beam element
%      sweepAngle       = structural sweep angle of element
%      aeroSweepAngle   = aerodynamic sweep angle of element
%      Omega            = rotor speed (Hz)
%      smOmega          = mode natural frequency for calculation of
%                         unsteady aero loads
%      airDensity       = freestream air density
%      aeroLoadsFlag    = flag to include or exclude aerodynamic loading
%                         0-exclude, 1-include
%      modalFlag        = flag to specify if element data is needed for
%                         modal or static analysis (1-modal, 0-static)
%      analysisType     = char for stability analysis type (F-standard
%                         rotating flutter analysis, P-parked flutter
%                         analysis, D-rotating divergence analysis)
%
%      output:
%      Ke               = element stiffness matrix
%      Me               = element mass matrix
%      Ce               = element damping matrix
%      Fe               = element load vector



%Set velocity and deactivate gyroscopic effects for parked analysis.
if(strcmpi(analysisType,'P'))
    Velocity = Omega;
    Omega = 0.0;
end


%calculate quad points
numGP = 4;
[xi,weight] = getGP(numGP);

%% Intitialize element matrices and vectors

%%Initialize element sub matrices and sub vectors
numNodesPerEl = length(x);
numDofPerNode = 6;

%stiffness sub-matrices
K11 = zeros(numNodesPerEl);
K12 = zeros(numNodesPerEl,numNodesPerEl*2);
K13 = zeros(numNodesPerEl,numNodesPerEl*2);
K14 = zeros(numNodesPerEl,numNodesPerEl);
K22 = zeros(numNodesPerEl*2,numNodesPerEl*2);
K23 = zeros(numNodesPerEl*2,numNodesPerEl*2);
K24 = zeros(numNodesPerEl*2,numNodesPerEl);
K33 = zeros(numNodesPerEl*2,numNodesPerEl*2);
K34 = zeros(numNodesPerEl*2,numNodesPerEl);
K44 = zeros(numNodesPerEl,numNodesPerEl);

%stress stiffening bending matrices
SS22 = K22;
SS33 = K33;

%spin softening sub-matrices
S11 = K11;
S12 = K12;
S13 = K13;
S14 = K14;
S22 = K22;
S23 = K23;
S24 = K24;
S33 = K33;
S34 = K34;
S44 = K44;

%load sub-vectors
F1 = zeros(numNodesPerEl,1);
F2 = zeros(numNodesPerEl*2,1);
F3 = F2;
F4 = F1;


if(modalFlag)
    %mass sub-matrices
    M11 = K11;
    M12 = K12;
    M13 = K13;
    M14 = K14;
    M22 = K22;
    M23 = K23;
    M24 = K24;
    M33 = K33;
    M34 = K34;
    M44 = K44;
    M42 = K24';
    M43 = M34';
    
    %coriolis sub-matrices
    C11 = K11;
    C12 = K12;
    C13 = K13;
    C14 = K14;
    C22 = K22;
    C23 = K23;
    C24 = K24;
    C32 = K23';
    C33 = K33;
    C34 = K34;
    C41 = K14';
    C42 = K24';
    C43 = K34';
    C44 = K44;
    
    %aerodynamic sub-matrices
    K22Aero = K22;
    K24Aero = K24;
    C24Aero = C42';
    C22Aero = C22;
    M22Aero = C22;
    M24Aero = K24;
    K44Aero = K44;
    K42Aero = K24';
    C42Aero = C42;
    C44Aero = K44;
    M42Aero = C42;
    M44Aero = K44;
end

%initialize overall element mass
Me = 0;
Ce = 0;

%%

%Convert frequencies from Hz to radians
smOmega = 2*pi*smOmega;
Omega = 2*pi*Omega/60.0;

%force rotor speed to unity for divergence analysis
if(strcmpi(analysisType,'D'))
    Omega = 1.0;
end

%Sort displacement vector
%Written for 2 node element with 6 dof per node

%calculate tranformation matrix for hub frame to element transformation
twistAvg = mean(sectionProps.twist);
lambda = calculateLambda(twistAvg*pi/180.0, coneAngle*pi/180.0, sweepAngle*pi/180.0);
lambdaSmall = lambda(1:3,1:3);

%convert element nodal displacements from hub frame from to element frame
disp = lambda*disp';

%sort element displacement vector into u, v, w, theta_x components
uNode = [disp(1) disp(7)];
vNode = [disp(2) disp(6) disp(8) disp(12)];
wNode = [disp(3) disp(5) disp(9) disp(11)];
theta_xNode = [disp(4)  disp(10)];

%Integration loop
for i=1:numGP
    %Calculate shape functions at quad point i
    [N,p_N_x,Jac] = calculateShapeFunctions(elementOrder,xi(i),xloc);
    [H,p_H_x,p_H_xx] = calculateHermiteShapeFunctions(xi(i),xloc,1.0);
    [H2,p_H_x2,p_H_xx2] = calculateHermiteShapeFunctions(xi(i),xloc,-1.0);
    integrationFactor = Jac * weight(i);
    
    %% Interpolate Element Properties at Quad Point
    
    EA   = interpolateVal(sectionProps.EA,N); %struct stiffness terms
    EIyy = interpolateVal(sectionProps.EIyy,N);
    EIzz = interpolateVal(sectionProps.EIzz,N);
    GJ   = interpolateVal(sectionProps.GJ,N);
    EIyz = interpolateVal(sectionProps.EIyz,N);
    Alpha = interpolateVal(sectionProps.Alpha,N);
    g1 = Alpha*sqrt(EIyy*GJ);
    alpha1   = interpolateVal(sectionProps.alpha1,N);
    alpha2   = interpolateVal(sectionProps.alpha2,N);
    alpha3   = interpolateVal(sectionProps.alpha3,N);
    alpha4   = interpolateVal(sectionProps.alpha4,N);
    alpha5   = interpolateVal(sectionProps.alpha5,N);
    alpha6   = interpolateVal(sectionProps.alpha6,N);
    
    rhoA   = interpolateVal(sectionProps.rhoA,N); %struct mass terms
    rhoIyy = interpolateVal(sectionProps.rhoIyy,N);
    rhoIzz = interpolateVal(sectionProps.rhoIzz,N);
    rhoIyz = interpolateVal(sectionProps.rhoIyz,N);
    
    a0gp = interpolateVal(sectionProps.a0,N);    %lift curve slope at gauss point
    Rgp  = interpolateVal(sectionProps.R,N) + hubRadius;     %rotor radius at gauss point
    bgp  = interpolateVal(sectionProps.b,N)/cos(sweepAngle*pi/180.0);     %semi chord at gauss point
    agp  = interpolateVal(sectionProps.a,N);     %elastic axis location at gauss point
    ycm = interpolateVal(sectionProps.ycm,N);
    zcm = interpolateVal(sectionProps.zcm,N);
    twistgp = interpolateVal(sectionProps.twist,N)*pi/180.0;
    acgp = interpolateVal(sectionProps.ac,N);
    
    rhoIyy = rhoIyy + rhoA*zcm^2;
    rhoIzz = rhoIzz + rhoA*ycm^2;
    rhoIyz = rhoIyz + rhoA*ycm*zcm;
    rhoJ = rhoIyy + rhoIzz;
    
    uprimegp = interpolateVal(uNode,p_N_x);
    vprimegp = interpolateVal(vNode,p_H_x);
    wprimegp = interpolateVal(wNode,p_H_x);
    
    xgp      = interpolateVal(x,N);
    ygp      = interpolateVal(y,N);
    zgp      = interpolateVal(z,N);
    
    axialForce = EA*uprimegp;
    %.... end interpolate value at quad points ........
    %%
    
    %calculate Uinfgp at quad point for rotating analysis
    if(strcmp(analysisType,'F') || strcmp(analysisType,'D'))
        Uinfgp = Rgp*Omega * cos(aeroSweepAngle*pi/180);
    end
    
    %calculate Uinfgp at quad point for parked analysis
    if(strcmp(analysisType,'P'))
        Uinfgp = Velocity * cos(aeroSweepAngle*pi/180);
    end
    
    %calculate reduced frequency and calculate Theodorsen function value
    kgp = smOmega*bgp/Uinfgp;
    Theogp = calculateTheo(kgp);
    if(strcmpi(analysisType,'D'))
        Theogp = 1.0;
    end
    
    %calculate moment arms for unsteady aerodynamics calculations
    d1gp = bgp*(agp+acgp);
    d2gp = bgp*(agp-0.5);
    
    %Calculate strutural stiffness sub matrices
    [K11] = calculateElType1(EA,integrationFactor,p_N_x,p_N_x,K11);
    [K12] = calculateElType2(0.5*EA*wprimegp*0,alpha3,integrationFactor,p_N_x,p_H_x,p_N_x,p_H_xx,K12);
    [K13] = calculateElType2(0.5*EA*vprimegp*0,alpha4,integrationFactor,p_N_x,p_H_x2,p_N_x,p_H_xx2,K13);
    [K14] = calculateElType1(alpha6,integrationFactor,p_N_x,p_N_x,K14);
    [K22] = calculateElType2(EIyy,0.5*EA*(wprimegp^2+vprimegp^2)*0,integrationFactor,p_H_xx,p_H_xx,p_H_x,p_H_x,K22);
    [K23] = calculateElType1(EIyz+alpha5,integrationFactor,p_H_xx2,p_H_xx2,K23);
    [K24] = calculateElType1(-g1,integrationFactor,p_H_xx,p_N_x,K24);
    [K33] = calculateElType2(EIzz,0.5*EA*(wprimegp^2+vprimegp^2)*0,integrationFactor,p_H_xx2,p_H_xx2,p_H_x2,p_H_x2,K33);
    [K34] = calculateElType1(alpha2,integrationFactor,p_H_xx2,p_N_x,K34);
    [K44] = calculateElType1(GJ,integrationFactor,p_N_x,p_N_x,K44);
    
    %Calculate stress stiffening bending matrices S22 and S33
    [SS22] = calculateElType1(axialForce,integrationFactor,p_H_x,p_H_x,SS22);
    [SS33] = calculateElType1(axialForce,integrationFactor,p_H_x2,p_H_x2,SS33);
    
    %% Calculate load vector
    %Calculate centrifugal load vector (only needed for static analysis)
    % define guass point load in hub coordinate system.
    
    %transform rotational velocity from hub to element frame
    OmegaVec = lambdaSmall*[0;0;Omega];
    O1 = OmegaVec(1);
    O2 = OmegaVec(2);
    O3 = OmegaVec(3);
    
    
    if(~modalFlag)
        %transform quad point position from hub to element frame
        posElCoordSys = lambdaSmall*[xgp;ygp;zgp];
        xl = posElCoordSys(1);
        yl = posElCoordSys(2);
        zl = posElCoordSys(3);
        
        %calculate quad point force/moments components
        f1 = rhoA*(xl*(O2^2+O3^2) - yl*O1*O2 - zl*O1*O3);
        f2 = rhoA*(zl*(O1^2+O2^2) - xl*O1*O3 - yl*O2*O3);
        f3 = rhoA*(yl*(O1^2+O3^2) - xl*O1*O2 - zl*O2*O3);
        f4 = rhoA*(xl*(zcm*O1*O2 - ycm*O1*O3) -yl*(ycm*O2*O3 +zcm*(O1^2+O3^2)) - zl*(zcm*O2*O3 +ycm*(O1^2+O2^2)));
        
        %calculate element load  sub-vectors
        [F1] = calculateVec1(f1,integrationFactor,N,F1);
        [F2] = calculateVec1(f2,integrationFactor,H,F2);
        [F3] = calculateVec1(f3,integrationFactor,H2,F3);
        [F4] = calculateVec1(f4,integrationFactor,N,F4);
    end
    
    %%
    %Spin softening/stiffening calculations
    [S11] = calculateElType1(-rhoA*(O2^2+O3^2),integrationFactor,N,N,S11);
    [S12] = calculateElType2(rhoA*O1*O3,rhoA*zcm*(O2^2+O3^2),integrationFactor,N,H,N,p_H_x,S12);
    [S13] = calculateElType2(rhoA*O1*O2,rhoA*ycm*(O2^2+O3^2),integrationFactor,N,H2,N,p_H_x2,S13);
    [S14] = calculateElType1(rhoA*(ycm*O1*O3-zcm*O1*O2),integrationFactor,N,N,S14);
    [S22] = calculateElType2(-rhoA*(O1^2+O2^2),-rhoIyy*(O2^2+O3^2),integrationFactor,H,H,p_H_x,p_H_x,S22);
    [S23] = calculateElType3(rhoA*O2*O3,-rhoIyz*(O1^2+O3^2),rhoA*(ycm*O1*O3-zcm*O1*O2),integrationFactor,H,H2,p_H_x,p_H_x2,p_H_x,H2,S23);
    [S24] = calculateElType2(rhoIyy*O1*O2-rhoIyz*O1*O3,-rhoA*(ycm*(O1^2+O2^2)+zcm*O2*O3),integrationFactor,p_H_x,N,H,N,S24);
    [S33] = calculateElType2(-rhoIzz*(O2^2+O3^2),-rhoA*(O1^2+O3^2),integrationFactor,p_H_x2,p_H_x2,H2,H2,S33);
    [S34] = calculateElType2(rhoIyz*O1*O2-rhoIzz*O1*O3,rhoA*(ycm*O2*O3 + zcm*(O1^2+O3^2)),integrationFactor,p_H_x2,N,H2,N,S34);
    [S44] = calculateElType1(-(rhoIyy*(O1^2+O3^2) +rhoIzz*(O1^2+O2^2) + 2*rhoIyz*O2*O3),integrationFactor,N,N,S44);
    
    if(modalFlag)
        %Calculate strutural mass sub matrices
        [M11] = calculateElType1(rhoA,integrationFactor,N,N,M11);
        [M12] = calculateElType1(rhoA*zcm,integrationFactor,N,p_H_x,M12);
        [M13] = calculateElType1(rhoA*ycm,integrationFactor,N,p_H_x2,M13);
        [M14] = calculateElType1(0,integrationFactor,N,N,M14);
        [M22] = calculateElType2(rhoIyy,rhoA,integrationFactor,p_H_x,p_H_x,H,H,M22);
        [M23] = calculateElType1(rhoIyz,integrationFactor,p_H_x,p_H_x2,M23);
        [M24] = calculateElType1(rhoA*ycm,integrationFactor,H,N,M24);
        [M33] = calculateElType2(rhoIzz,rhoA,integrationFactor,p_H_x2,p_H_x2,H2,H2,M33);
        [M34] = calculateElType1(rhoA*zcm,integrationFactor,H2,N,M34);
        [M44] = calculateElType1(rhoJ,integrationFactor,N,N,M44);
        M42=M24';
        
        
        
        %damping matrix calculations (coriolis)
        C12 = calculateElType1(2*rhoA*O2,integrationFactor,N,H,C12);
        C13 = calculateElType1(-2*rhoA*O3,integrationFactor,N,H2,C13);
        C14 = calculateElType1(2*rhoA*(ycm*O2+zcm*O3),integrationFactor,N,N,C14);
        C22 = calculateElType2(2*rhoA*zcm*O2,-2*rhoA*zcm*O2,integrationFactor,H,p_H_x,p_H_x,H,C22);
        C23 = calculateElType2(2*rhoA*O1,2*rhoA*(zcm*O3-ycm*O2),integrationFactor,H,H2,p_H_x,H2,C23);
        C24 = calculateElType2(-(2*rhoIyy*O3+2*rhoIyz*O2),-2*rhoA*zcm*O1,integrationFactor,p_H_x,N,H,N,C24);
        C33 = calculateElType2(-2*rhoA*ycm*O3,2*rhoA*ycm*O3,integrationFactor,p_H_x2,H2,H2,p_H_x2,C33);
        C34 = calculateElType2(-(2*rhoIzz*O2 + 2*rhoIyz*O3),-2*rhoA*ycm*O1,integrationFactor,p_H_x2,N,H2,N,C34);
        
        %enforce skew symmetry of coriolis matrix
        C21 = -C12';
        C31 = -C13';
        C41 = -C14';
        C32 = -C23';
        C42 = -C24';
        C43 = -C34';
        
        %% ...calculate aerodynamic matrices........
        if(aeroLoadsFlag==1 || aeroLoadsFlag==3)
            [K24Aero] = calculateElType1(-a0gp*airDensity*Uinfgp^2*bgp*Theogp,integrationFactor,H,N,K24Aero);
            [C24Aero] = calculateElType1(-0.5*a0gp*airDensity*Uinfgp*bgp^2*(1+Theogp*(1-2*agp)),integrationFactor,H,N,C24Aero);
            [C22Aero] = calculateElType1(a0gp*airDensity*Uinfgp*bgp*Theogp,integrationFactor,H,H,C22Aero);
            [M22Aero] = calculateElType1(0.5*a0gp*airDensity*bgp^2,integrationFactor,H,H,M22Aero);
            [M24Aero] = calculateElType1(0.5*a0gp*airDensity*bgp^3*agp,integrationFactor,H,N,M24Aero);
            
            [K44Aero] = calculateElType1(-a0gp*airDensity*Uinfgp^2*bgp*d1gp*Theogp,integrationFactor,N,N,K44Aero);
            [C42Aero] = calculateElType1(a0gp*airDensity*d1gp*Uinfgp*bgp*Theogp,integrationFactor,N,H,C42Aero);
            [C44Aero] = calculateElType1(-0.5*a0gp*airDensity*Uinfgp*bgp^2*(d1gp*(Theogp*(1-2.0*agp))+d2gp),integrationFactor,N,N,C44Aero);
            [M42Aero] = calculateElType1(0.5*a0gp*airDensity*bgp^3*agp,integrationFactor,N,H,M42Aero);
            [M44Aero] = calculateElType1(0.5*a0gp*airDensity*bgp^4*(0.125+agp^2),integrationFactor,N,N,M44Aero);
            
            % Chris Kelley (implement Greenberg theory for surge (edgewise)
            % motion. Torque 2020
            if(aeroLoadsFlag==3)
                alpha = 6/180*pi; % Depends on steady angle of attack
                sigma = 1; % amplitude ratio y_dot/W 
                [M21Aero] = calculateElType1(0.5*a0gp*airDensity*bgp^2*alpha*sigma,integrationFactor,H,H,M21Aero)
                [M41Aero] = calculateElType1(-0.5*a0gp*airDensity*bgp^3*alpha*agp*sigma,integrationFactor,N,H,M41Aero)
                [C21Aero] = calculateElType1(a0gp*airDensity*bgp*Uinfgp*alpha*(1+Theogp)*sigma,integrationFactor,H,H,C21Aero)
                [C41Aero] = calculateElType1(-a0gp*airDensity*bgp^2*alpha*(agp+0.5)*Uinfgp*alpha*(1+Theogp)*sigma,integrationFactor,N,H,C41Aero)
                
            end
            
            
            
        end
        
        if(aeroLoadsFlag==2)
            %This is a real valued aeroelastic representation from
            %Wright and Cooper
            %This version assumes aerodynamic center at quarter chord of
            %airfoil
            Fgp = real(Theogp);
            Ggp = imag(Theogp);
            lcsrat = a0gp/(2*pi);
            Fgp = Fgp*lcsrat;
            Ggp = Ggp*lcsrat;
            
            if(Uinfgp==0)
                kgp = 1;
            end
            
            lz = -2*pi*(-0.5*kgp^2-Ggp*kgp);
            lzdot = -2*pi*Fgp;
            ltheta = 2*pi*(0.5*kgp^2*agp + Fgp - Ggp*kgp*(0.5-agp));
            lthetadot = 2*pi*(.5 + Fgp*(.5-agp) + Ggp/kgp);
            
            mz = -2*pi*(-0.5*kgp^2*agp-kgp*(agp+.5)*Ggp);
            mzdot = -2*pi*(agp+0.5)*Fgp;
            mtheta = 2*pi*(0.5*kgp^2*(1.0/8.0+agp^2)+Fgp*(agp+0.5)-kgp*Ggp*(agp+0.5)*(0.5-agp));
            mthetadot = 2*pi*(-0.5*kgp*(0.5-agp) + kgp*Fgp*(agp+0.5)*(.5-agp)+Ggp/kgp*(agp+0.5));
            
            
            
            k22fac = airDensity*Uinfgp^2*lz;
            k24fac = airDensity*Uinfgp^2*bgp*ltheta;
            k42fac = airDensity*Uinfgp^2*bgp*mz;
            k44fac = airDensity*Uinfgp^2*bgp^2*mtheta;
            
            c22fac = airDensity*Uinfgp*bgp*lzdot;
            c24fac = airDensity*Uinfgp*bgp^2*lthetadot;
            c42fac = airDensity*Uinfgp*bgp^2*mzdot;
            c44fac = airDensity*Uinfgp*bgp^3*mthetadot;
            
            
            
            [K22Aero] = calculateElType1(-k22fac,integrationFactor,H,H,K22Aero);
            [K24Aero] = calculateElType1(-k24fac,integrationFactor,H,N,K24Aero);
            [C24Aero] = calculateElType1(-c24fac,integrationFactor,H,N,C24Aero);
            [C22Aero] = calculateElType1(-c22fac,integrationFactor,H,H,C22Aero);
            
            [K42Aero] = calculateElType1(-k42fac,integrationFactor,N,H,K42Aero);
            [K44Aero] = calculateElType1(-k44fac,integrationFactor,N,N,K44Aero);
            [C42Aero] = calculateElType1(-c42fac,integrationFactor,N,H,C42Aero);
            [C44Aero] = calculateElType1(-c44fac,integrationFactor,N,N,C44Aero);
            
        end
        
        
        %........................................
    end
end

%combine structural stiffness and spin softening matrices for
%rotational analysis
if(~strcmpi(analysisType,'D'))
    K11 = K11 + S11;
    K21 = 2*K12' + S12';
    K12 = K12 + S12;
    K31 = 2*K13' + S13';
    K13 = K13 + S13;
    K14 = K14 + S14;
    K22 = K22 + S22 + SS22;
    K23 = K23 + S23;
    K24 = K24 + S24;
    K33 = K33 + S33 + SS33;
    K34 = K34 + S34;
    K44 = K44 + S44;
end
%---------------------------------------------

%aerodynamic contributions--------------------
K42 = K24';

%combine aerodynamic matrices with overall element matrices
if(aeroLoadsFlag && ~strcmpi(analysisType,'D'))
    K22 = K22 + K22Aero; %used for real valued rep
    K24 = K24 + K24Aero; %K is symmmetric before adding aero terms
    K42 = K42 + K42Aero;
    K44 = K44 + K44Aero;
    
    C22 = C22 + C22Aero; %C is non-sym before adding aero terms
    C24 = C24 + C24Aero;
    C42 = C42 + C42Aero;
    C44 = C44Aero;
    if(aeroLoadsFlag==3)
    C21 = C21 + C21Aero;
    C41 = C41 + C41Aero;
    end
    
    M22 = M22 + M22Aero; %M is symmetric before adding aero terms
    M43 = M34';          %MijAero not used in real valued rep
    M24 = M24 + M24Aero;
    M42 = M24'+ M42Aero;
    M44 = M44 + M44Aero;
    if(aeroLoadsFlag==3)
    M21 = M21 + M21Aero;
    M41 = M41 + M41Aero;
    end
    
end


%create overall system matrices for diveregence analysis
if(strcmpi(analysisType,'D'))
    %This code assumes aeroLoadsFlags = 1
    %redefine K to have pure structural stiffness
    K21 = 2*K12';
    K31 = 2*K13';
    K42 = K24';
    
    
    %overwrite mass matrix with Omega^2Coefficient matrix
    M11 = S11;
    M21 = S12';
    M12 = S12;
    M31 = S13';
    M13 = S13;
    M14 = S14;
    M22 = S22;
    M23 = S23;
    M24 = S24 + K24Aero;
    M42 = S24';
    M33 = S33;
    M34 = S34;
    M44 = S44 + K44Aero;
end
%---------------------------------------------

%map element matrices for conventional DOF numbering
[Ke] = mapMatrixNonSym(K11,K12,K13,K14,...
    K21,K22,K23,K24,...
    K31,K23',K33,K34,...
    K14',K42,K34',K44);

if(modalFlag)
    [Ce] = mapMatrixNonSym(C11,C12,C13,C14,...
        C21,C22,C23,C24,...
        C31,C32,C33,C34,...
        C41,C42,C43,C44);
    
    if(strcmpi(analysisType,'D'))
        [Me] = -1.0*mapMatrixNonSym(M11,M12,M13,M14,...
            M21,M22,M23,M24,...
            M31,M23',M33,M34,...
            M14',M42,M34',M44);
    else
        [Me] = mapMatrixNonSym(M11,M12,M13,M14,...
            M12',M22,M23,M24,...
            M13',M23',M33,M34,...
            M14',M42,M43,M44);
    end
    
end

%map element vector for conventional DOF numbering
[Fe] = mapVector(F1,F2,F3,F4);

% transform matrices and load vector from element to hub frame
lambdaTran = lambda';

if(modalFlag)
    Me = lambdaTran*Me*lambda;
    Ce = lambdaTran*Ce*lambda;
end

Ke = lambdaTran*Ke*lambda;
Fe = lambdaTran*Fe;


end

function [valGP] = interpolateVal(valNode,N)
% This function interpolates values of an element at the gauss point
%
%    input:
%    valNode    = nodal values
%    N          = shape function evaluated at gausspoint
%
%   output:
%   valGP       = value at gauss point

valGP = 0.0;
for i=1:length(N)
    valGP = valGP + N(i)*valNode(i);
end
end

function [K] = calculateElType1(fac,integrationFactor,N1,N2,K)
% This function performs a generic element matrix calculation
%
%    input:
%    fac               = integrand shape function coefficient
%    integrationFactor = element length associated with gauss point
%    N1,N2             = generic shape function values
%    K                 = element matrix to be added to
%
%    output:
%    K                 = output element matrix

K = K + N1'*N2*fac*integrationFactor;

end

function [K] = calculateElType2(fac1,fac2,integrationFactor,N1,N2,N3,N4,K)
% This function performs a generic element matrix calculation
%
%    input:
%    fac1, fac2        = integrand shape function coefficients
%    integrationFactor = element length associated with gauss point
%    N1,N2,N3,N4       = generic shape function values
%    K                 = element matrix to be added to
%
%    output:
%    K                 = output element matrix

%N1 and N3 must be same length
%N2 and N4 must be same length

K = K + (N1'*N2*fac1 + N3'*N4*fac2)*integrationFactor;
end

function [K] = calculateElType3(fac1,fac2,fac3,integrationFactor,N1,N2,N3,N4,N5,N6,K)
% This function performs a generic element matrix calculation
%
%    input:
%    fac1, fac2, fac3  = integrand shape function coefficients
%    integrationFactor = element length associated with gauss point
%    N1,N2,N3,N4,N5,N6 = generic shape function values
%    K                 = element matrix to be added to
%
%    output:
%    K                 = output element matrix

%N1 and N3 and N5 must be same length
%N2 and N4 and N6 must be same length

K = K + (fac1*N1'*N2 + fac2*N3'*N4 + fac3*N5'*N6)*integrationFactor;
end

function [F] = calculateVec1(f,integrationFactor,N,F)
% This function performs a generic element vector calculation
%
%    input:
%    f                 = integrand shape function coefficients
%    integrationFactor = element length associated with gauss point
%    N                 = generic shape function values
%    F                 = element vector to be added to
%
%   output:
%   F                  = output element vector
F = F + f*N'*integrationFactor;

end

function [Kel] = mapMatrix(K11,K12,K13,K14,K22,K23,K24,K33,K34,K44)
%This function maps submatrices to entries in a symmetric matrix
%
%    input:
%    Kij   = input submatrices
%
%    output:
%    Kel    = populated full element matrix with mapped submatrix values

%create a full stiffness matrix with DOF ordering as derived from equations
Ktemp = [K11, K12,K13,K14;
    K12',K22,K23,K24;
    K13',K23',K33,K34;
    K14',K24',K34',K44];
[a]=length(Ktemp);
Kel = zeros(a);

%declare map
map = [1, 7, 3, 5, 9, 11, 2, 6, 8, 12, 4, 10];



% Kel = zeros(a,a);

%map to FEA numbering
for i=1:a
    I=map(i);
    for j=1:a
        J=map(j);
        Kel(I,J) = Ktemp(i,j);
    end
end

end


function [Kel] = mapMatrixNonSym(K11,K12,K13,K14,K21,K22,K23,K24,...
    K31,K32,K33,K34,...
    K41,K42,K43,K44)
%This function maps submatrices to entries in non-symmetric matrix
%
%    input:
%    Kij   = input submatrices
%
%    output:
%    Kel    = populated full element matrix with mapped submatrix values

%create a full stiffness matrix with DOF ordering as derived from equations
Ktemp = [K11, K12,K13,K14;
    K21,K22,K23,K24;
    K31,K32,K33,K34;
    K41,K42,K43,K44];

a=length(Ktemp);
Kel = zeros(a);

%declare map
map = [1, 7, 3, 5, 9, 11,...
    2, 6, 8, 12, 4, 10];

%map to FEA numbering
for i=1:a
    I=map(i);
    for j=1:a
        J=map(j);
        Kel(I,J) = Ktemp(i,j);
    end
end

end


function [Fel] = mapVector(F1,F2,F3,F4)

Ftemp = [F1;F2;F3;F4];
a=length(Ftemp);
Fel=zeros(a,1);

%declare map
map = [1, 7, 3, 5, 9, 11,...
    2, 6, 8, 12, 4, 10];

%map to FEA numbering
for i=1:a
    I=map(i);
    Fel(I) = Ftemp(i);
end

end
