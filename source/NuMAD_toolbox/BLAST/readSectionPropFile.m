function [model] = readSectionPropFile(bladefn,adfn,bladeLength,model)
%readSectionPropsFile Reads blade section properties from file
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [model] = readSectionPropFile(bladefn,adfn,bladeLength,model)
%                    
%   This function opens and reads the blade section property file and
%   aerodyn file for reading geometric, structural, and aerodynamic
%   properties of the blade. Differences between the FAST and BLAST
%   coordinate systems are accounted for in this function through this
%   transformation.
%
%      input:
%      bladefn      = blade filename
%      adfn         = AeroDyn filename
%      bladeLength  = length of blade
%      model        = model struct to have section properties appended
%
%      output:
%      model        = model struct output with section property information

    blade = readFastBlade(bladefn);
    ad = readFastAD(adfn);  
    model.airDensity = ad.Rho;
    numSections = length(blade.prop.EAStff);
    
    model.numElements = numSections-1;
    
    R =  blade.prop.BlFract.*bladeLength;
    model.R=R;
    aeroR = ad.RNodes-model.hubRadius;
    EA = blade.prop.EAStff;
    EIyy = blade.prop.FlpStff;
    EIzz = blade.prop.EdgStff;
    GJ = blade.prop.GJStff;
    EIyz = zeros(numSections,1);
    Alpha = blade.prop.Alpha;
    alpha1 = zeros(numSections,1);
    alpha2 = zeros(numSections,1);
    alpha3 = zeros(numSections,1);
    alpha4 = zeros(numSections,1);
    alpha5 = zeros(numSections,1);
    alpha6 = zeros(numSections,1);
    rhoA = blade.prop.BMassDen;
    rhoIyy = blade.prop.FlpIner;
    rhoIzz = blade.prop.EdgIner;
    rhoJ = rhoIyy + rhoIzz;
    rhoIyz = zeros(numSections,1);
    structuralTwist = blade.prop.StrcTwst;
    bAeroNode = ad.Chord.*0.5;
    a = blade.prop.EdgEAOf;
    ycm = blade.prop.EdgcgOf;
    zcm = blade.prop.FlpcgOf;
    airfoilNum = ad.NFoil;
    preSweep = blade.prop.PreswpRef;
    preCurve = blade.prop.PrecrvRef;
    
    fstac = blade.prop.AeroCent;
    
    model.preSweep = -preSweep;
    model.edgeEAOffset = -a;
    model.preCurve = preCurve;
    model.flapEAOffset = blade.prop.FlpEAOf;
    
    model.mesh.y = -blade.prop.EdgEAOf;
    model.mesh.z = blade.prop.FlpEAOf;
      
    LCS = model.LCS;
    for i=1:length(aeroR)
        a0AeroNode(i) = LCS(airfoilNum(i));
    end
    
    [b,a0] = interpolateAeroDataAtBladeSections(R,aeroR,bAeroNode,a0AeroNode);
               
    ptchAx=model.ptchAx;
    ac = (0.5-(fstac - 0.25 + ptchAx))*2.0; %converts from FAST AC to actual AC and expresses in semi-chord fraction

%     ycm = 2.*ycm - 1.0;
     ycm = -(ycm-a);
     zcm = (zcm-blade.prop.FlpEAOf);
    aDistAftLE = (a + ptchAx.*2.*b)./(2.*b); %chord fractions aft of LE
     %convert to semi-chord fraction from half-chord (positive aft - towards TE)
    a  = 2.*aDistAftLE - 1.0;
    
    [model.sweepAngle,model.coneAngle] = calculateSweepAngle(-blade.prop.EdgEAOf,blade.prop.FlpEAOf,R,-preSweep,preCurve);
    model.aeroSweepAngle = calculateAeroSweepAngle(R,-preSweep);
    ptchAx = ptchAx - preSweep; 
    plotBladePlanform(R,ptchAx,aDistAftLE,b,ycm,ac);


    
%     [len,dum] = size(data);
    
    for i=1:numSections-1
        sectionPropsArray{i}.R = [R(i), R(i+1)];
        
        sectionPropsArray{i}.EA = [EA(i), EA(i+1)];
        sectionPropsArray{i}.EIyy = [EIyy(i), EIyy(i+1)];
        sectionPropsArray{i}.EIzz = [EIzz(i), EIzz(i+1)];
        sectionPropsArray{i}.GJ = [GJ(i), GJ(i+1)];
        sectionPropsArray{i}.EIyz = [EIyz(i), EIyz(i+1)];
        sectionPropsArray{i}.Alpha = [Alpha(i), Alpha(i+1)];
        sectionPropsArray{i}.alpha1 = [alpha1(i), alpha1(i+1)];
        sectionPropsArray{i}.alpha2 = [alpha2(i), alpha2(i+1)];
        sectionPropsArray{i}.alpha3 = [alpha3(i), alpha3(i+1)];
        sectionPropsArray{i}.alpha4 = [alpha4(i), alpha4(i+1)];
        sectionPropsArray{i}.alpha5 = [alpha5(i), alpha5(i+1)];
        sectionPropsArray{i}.alpha6 = [alpha6(i), alpha6(i+1)];
        
        sectionPropsArray{i}.rhoA = [rhoA(i), rhoA(i+1)];
        sectionPropsArray{i}.rhoIyy = [rhoIyy(i), rhoIyy(i+1)];
        sectionPropsArray{i}.rhoIzz = [rhoIzz(i), rhoIzz(i+1)];
        sectionPropsArray{i}.rhoJ = [rhoJ(i), rhoJ(i+1)];
        sectionPropsArray{i}.rhoIyz = [rhoIyz(i), rhoIyz(i+1)];
        
        sectionPropsArray{i}.b = [b(i), b(i+1)];
        sectionPropsArray{i}.a = [a(i), a(i+1)];
        sectionPropsArray{i}.ac =[ac(i), ac(i+1)];
        sectionPropsArray{i}.ycm = [ycm(i), ycm(i+1)];
        sectionPropsArray{i}.zcm = [zcm(i), zcm(i+1)];
        sectionPropsArray{i}.twist=[structuralTwist(i), structuralTwist(i+1)];
        sectionPropsArray{i}.a0 = [a0(i), a0(i+1)];
    end
   
model.sectionPropsArray = sectionPropsArray;

end

function [b,a0] = interpolateAeroDataAtBladeSections(R,aeroR,bAeroNode,a0AeroNode)
% This function interpolates data from aero domain to structural domain
%
%    input:
%    R          = spanwise nodes in structural domain
%    aeroR      = spanwise nodes in aero domain
%    bAeroNode  = semi chord at aero nodes
%    a0AeroNode = lift curve slope at aero nodes
%
%    output:
%    b          = semi chord at structural nodes
%    a0         = lift curve slope at structural nodes

    %interpolate from aeronodes section
          %set root point and tip point
     if(aeroR(1) ~= 0.0)
         aeroR = [R(1); aeroR];
         bAeroNode = [bAeroNode(1); bAeroNode];
         a0AeroNode = [a0AeroNode(1) a0AeroNode];
     end
     
     lenAeroR=length(aeroR);
     if(aeroR(lenAeroR) < R(length(R)))
         
         bslope = (bAeroNode(lenAeroR)-bAeroNode(lenAeroR-1))/(aeroR(lenAeroR)-aeroR(lenAeroR-1));
         a0slope = (a0AeroNode(lenAeroR)-a0AeroNode(lenAeroR-1))/(aeroR(lenAeroR)-aeroR(lenAeroR-1));
         
         aeroR = [aeroR; R(length(R))];
         delRTip = aeroR(lenAeroR+1)-aeroR(lenAeroR);
         bAeroTip = bAeroNode(lenAeroR) + bslope*delRTip;
         a0AeroTip = a0AeroNode(lenAeroR) + a0slope*delRTip;
         
         
         bAeroNode = [bAeroNode; bAeroTip];
         a0AeroNode = [a0AeroNode a0AeroTip];
     end
     
     b = interp1(aeroR,bAeroNode,R);
     a0 = interp1(aeroR,a0AeroNode,R);
     % end interpolate from aero nodes section
    

end


function [lambda,cone] = calculateSweepAngle(elasticAxisEdge,elasticAxisFlap,sectionLoc,preSweep,preCurve)
% This function calculates the structural sweep angle
%
%    input:
%    elasticAxis = location of elastic axis at nodes
%    sectionLoc  = spanwise section locations
%    semiChord   = semi chord values at nodes
%    preSweep    = pre sweep value at nodes
%
%    output:
%    lambda      = structural sweep angle (degrees)
elasticAxisEdge = elasticAxisEdge + preSweep;
elasticAxisFlap  = elasticAxisFlap + preCurve;
len = length(elasticAxisEdge);
lambda = zeros(len-1,1);
cone = zeros(len-1,1);

for i=1:len-1
    delElastAxEdge(i) = elasticAxisEdge(i+1) - elasticAxisEdge(i);
    delElastAxFlap(i) = elasticAxisFlap(i+1) - elasticAxisFlap(i);

    delSectionLoc(i) = sectionLoc(i+1)-sectionLoc(i);
    
    
    lambda(i) = asin(delElastAxEdge(i)/sqrt(delSectionLoc(i)^2 + delElastAxEdge(i)^2));
    cone(i)   = -asin(delElastAxFlap(i)/sqrt(delSectionLoc(i)^2 + delElastAxEdge(i)^2+delElastAxFlap(i)^2));
    lambda(i) = lambda(i)*180.0/pi;
    cone(i) = cone(i)*180.0/pi;
end

end

function [aeroLambda] = calculateAeroSweepAngle(sectionLoc,preSweep)
% This function calculates the structural sweep angle
%
%    input:
%    sectionLoc  = spanwise section locations
%    preSweep    = pre sweep value at nodes
%
%    output:
%    lambda      = aero sweep angle (degrees)

len = length(sectionLoc);
aeroLambda = zeros(len-1,1);

for i=1:len-1
    delpreSwp = preSweep(i+1) - preSweep(i);
    delSectionLoc = sectionLoc(i+1)-sectionLoc(i);
    
    aeroLambda(i) = asin(delpreSwp/delSectionLoc);
    aeroLambda(i) = aeroLambda(i)*180.0/pi;
end

end

function plotBladePlanform(R,ptchAx,aDistAftLE,b,ycm,ac)
% This function plots the blade planform and axes
%
%    input:
%    R          = spanwise location of blade sections (nodes)
%    ptchAx     = pitch axis values at nodes
%    aDistAFTLE = distance elastic axis is aft of leading edge at nodes
%    b          = semi chord value at nodes
%    ycm        = edgewise mass center offset from elastic axis at nodes
%    ac         = aerodynamic center value at nodes
%
%    output:
%    none - No explicit output from function, plot is generated
    figure(100);
    
    plot(R,-(ptchAx-aDistAftLE).*2.*b,'--r',... %Elastic Axis
         R,-((ptchAx-aDistAftLE).*2.*b+ycm),'-go',... %Mass Center
         R,-(2.*ptchAx-1).*b,'--k',... %Semi-chord
         R,-((2.*ptchAx-(1.0-ac)).*b),'-cd',...%Aero Center
         R,-ptchAx.*2.*b,'k',... %LE
         R,-(ptchAx-1).*2.*b,'k','LineWidth',2.0); %TE
        
     
    legend('Elastic Axis','Mass Center',...
            'Semi Chord','Aero Center');

    grid on
    axis equal
    hold off
    
    %Label LE and TE on plot instead of in legend
    mTextBox = uicontrol('style','text');
    set(mTextBox,'String','Leading Edge','Position',[80 160 110.0 20.0],'FontSize',12);
    mTextBox2 = uicontrol('style','text');
    set(mTextBox2,'String','Trailing Edge','Position',[80 250 120.0 20.0],'FontSize',12);
    xlabel('Span Location (m)');
    ylabel('Edgewise Axes Location (m)');
end

