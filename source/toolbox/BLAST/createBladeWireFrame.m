function [B1,B2,B3l,B3u,eaLoc2,eaLoc3] = createBladeWireFrame(bladefn,adfn,bladeLength,Psi,Theta,SweepAng,offset,plotFlag)
%createBladeWireFrame Creates wireframe data and 3D plot of blade
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [B1,B2,B3l,B3u,eaLoc2,eaLoc3] = createBladeWireFrame(bladefn,adfn,
%                                   bladeLength,Psi,Theta,SweepAng,offset,
%                                   plotFlag)
%                    
%   This function accepts blade and aerodyn data files, blade parameters,
%   blade orientation and offset. Data for the blade wireframe is output
%   and the blade wireframe is plotted.
%
%      input:
%      bladefn      = blade data file filename
%      adfn         = aerodyn data file filename
%      bladeLength  = length of blade
%      Psi          = azimuth angle (degrees) about the h3 axis
%      Theta        = cone angle (degrees) for blade orientation
%      SweepAng     = local sweep angle (degrees) for blade orientation
%      offset       = 3x1 vector for translating component from origin
%      plotFlag     = flag to control plotting: 1-plots, 0-no plots
%
%      output:
%      B1           = h1 coordinates of blade wireframe
%      B2           = h2 coordinates of blade wireframe
%      B3l          = h3 coordinates of lower wireframe
%      B3u          = h3 coordinates of upper wireframe
%      eaLoc2       = h2 coordinates of elastic axis for blade
%      eaLoc3       = h3 coordinates of elastic axis for blade

    blade = readFastBlade(bladefn);
    ad = readFastAD(adfn);  
    
    stationLoc = blade.prop.BlFract .* bladeLength;
    stationTwist = blade.prop.StrcTwst.*pi/180;
    
    [b,nfoil] = interpolateAeroDataAtBladeSections(stationLoc,ad.RNodes,0.5.*ad.Chord,ad.NFoil);
    
    deltaTheta = 15*pi/180;
    theta=[0:deltaTheta:pi];
    leny=length(theta);
    lenb1 = length(stationLoc);
        
    for i=1:lenb1
        if(nfoil(i) == 1)
            thicknessRatio = 1.0;
        elseif(nfoil(i) == 2)
            thicknessRatio = 0.6;
        elseif(nfoil(i) == 3)
            thicknessRatio = 0.3;
        else
            thicknessRatio = 0.2;
        end
        
        B1(i,:) = stationLoc(i)*ones(1,leny);
        B2(i,:) = b(i).*cos(theta) - blade.prop.PreswpRef(i);
        B3u(i,:) = b(i).*sin(theta).*thicknessRatio + blade.prop.PrecrvRef(i);
        B3l(i,:) = -b(i).*sin(theta).*thicknessRatio + blade.prop.PrecrvRef(i);
        
        eaLoc2(i) = blade.prop.PreswpRef(i) + blade.prop.EdgEAOf(i);
        eaLoc3(i) = blade.prop.PrecrvRef(i) + blade.prop.FlpEAOf(i);
        [~,B3u(i,:)] = twistSection(stationTwist(i),B2(i,:),B3u(i,:));
        [B2(i,:),B3l(i,:)] = twistSection(stationTwist(i),B2(i,:),B3l(i,:));
    end
    
    [H1u,H2u,H3u] = rigidRotationBlade(B1,B2,B3u,Psi,Theta,SweepAng);
    [H1l,H2l,H3l] = rigidRotationBlade(B1,B2,B3l,Psi,Theta,SweepAng);
    H1l=flipud(H1l);
    H2l=flipud(H2l);
    H3l=flipud(H3l);
    
    H1 = cat(1,H1u,H1l);
    H2 = cat(1,H2u,H2l);
    H3 = cat(1,H3u,H3l);
    
    H1 = H1 + offset(1);
    H2 = H2 + offset(2);
    H3 = H3 + offset(3);
    
   if(plotFlag)     
    hold on;
%     axis equal;
    mesh(B1,B2,B3u,'EdgeColor','black');
    mesh(B1,B2,B3l,'EdgeColor','black');
%     mesh(H1,H2,H3,'EdgeColor','red');
%     view(3);
    alpha(0.3);
   end

end

function [B2new,B3new] = twistSection(twistAngle,B2,B3)
% This function twists a blade section
%
%    input:
%    twistAngle    = angle of twist
%    B2            = 2 coordinates of section
%    B3            = 3 coordinates of section
%
%    output:
%    B2new         = new 2 coordinates of section
%    B3new         = new 3 coordinates of section

    ct = cos(twistAngle);
    st = sin(twistAngle);
    C = [1 0 0;0 ct st;0 -st ct]';
    
    B2new = C(2,2).*B2 + C(2,3).*B3;
    B3new = C(3,2).*B2 + C(3,3).*B3;
    
end

function [H1,H2,H3] = rigidRotationBlade(B1,B2,B3,Psi,Theta,SweepAng)
%Performs a rigid rotation of the blade through a 3-2-3 Euler sequence
%
%    input:
%    B1        = 1 coordinates of body in body frame
%    B2        = 2 coordinates of body in body frame
%    B3        = 3 coordinates of body in body frame
%    Psi       = 1st 3 rotation angle
%    Theta     = 2 rotation angle
%    SweepAng  = 2nd 3 rotation angle
%
%    output:
%    H1        = 1 coordiantes of oriented body in hub frame
%    H2        = 2 coordinates of oriented body in hub frame
%    H3        = 3 coordinates of oriented body in hub frame

        %peform rotation about h3 axis
        Psi = Psi*pi/180.0;
        Theta = Theta*pi/180.0;
        SweepAng = SweepAng*pi/180.0;
        
        M3 = [cos(Psi) sin(Psi) 0; -sin(Psi) cos(Psi) 0; 0 0 1];
        M2 = [cos(Theta) 0 -sin(Theta);0 1 0; sin(Theta),0,cos(Theta)];
        M3f = [cos(SweepAng) sin(SweepAng) 0; -sin(SweepAng) cos(SweepAng) 0; 0 0 1];
        
        C=(M3f*M2*M3)';
    
        H1 = C(1,1).*B1 + C(1,2).* B2 + C(1,3).*B3;
        H2 = C(2,1).*B1 + C(2,2).* B2 + C(2,3).*B3;
        H3 = C(3,1).*B1 + C(3,2).* B2 + C(3,3).*B3;
end

function [b,nfoil] = interpolateAeroDataAtBladeSections(R,aeroR,bAeroNode,nFoilAeroNode)
% This function interpolates data from the aero domain to the structural
% blade domain
%
%  input:
%  R              = blade section locations for structural domain
%  aeroR          = blade section locations for aero domain
%  bAeroNode      = semi chord at aero nodes
%  nFoilAeroNodes = airfoil number at aero node
%
%  output:
%  b              = semi chord projected on structural domain
%  nfoil          = air foil number projected on structural domain

    %interpolate from aeronodes section
          %set root point and tip point
     if(aeroR(1) ~= 0.0)
         aeroR = [R(1); aeroR];
         bAeroNode = [bAeroNode(1); bAeroNode];
     end
     
     lenAeroR=length(aeroR);
     if(aeroR(lenAeroR) < R(length(R)))
         
         bslope = (bAeroNode(lenAeroR)-bAeroNode(lenAeroR-1))/(aeroR(lenAeroR)-aeroR(lenAeroR-1));
                  
         aeroR = [aeroR; R(length(R))];
         delRTip = aeroR(lenAeroR+1)-aeroR(lenAeroR);
         bAeroTip = bAeroNode(lenAeroR) + bslope*delRTip;
         
         bAeroNode = [bAeroNode; bAeroTip];
         nFoilAeroNode = [nFoilAeroNode(1);nFoilAeroNode;nFoilAeroNode(length(nFoilAeroNode))];
     end
     
     b = interp1(aeroR,bAeroNode,R);
     nfoil = interp1(aeroR,nFoilAeroNode,R,'nearest');
     % end interpolate from aero nodes section
 
end

