classdef shellRegion
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type = '';
        keyPts = [];
        edgeEls = [];
    end
    
    methods
        function obj = shellRegion(regType,keyPoints,numEdgeEls)
            obj.type = regType;
            obj.keyPts = keyPoints;
            obj.edgeEls = numEdgeEls;
        end
        
        function [nodes,elements] = createShellMesh(obj,elType,method)
            if(contains(method,'structured'))
                if(contains(obj.type,'quad') || contains(obj.type,'cyl'))
                    xNodes = max([obj.edgeEls(1),obj.edgeEls(3)],[],'all') + 1;
                    yNodes = max([obj.edgeEls(2),obj.edgeEls(4)],[],'all') + 1;
                    boundaryNodes = [linspace(-1,1,xNodes)',-1*ones(xNodes,1)];
                    boundaryEdges = [1:(xNodes-1);2:xNodes]';
                    mesh = NuMesh2D(boundaryNodes,boundaryEdges);
                    [nodes,elements] = mesh.createSweptMesh('in_direction',[0,1],2,(yNodes-1),0,[]);
                    ndElim = zeros(length(nodes),1);
                    if(obj.edgeEls(1) < obj.edgeEls(3))
                        nds2Elim = obj.edgeEls(3) - obj.edgeEls(1);
                        rowstep = ceil((xNodes-nds2Elim)/(nds2Elim+1)) + 1;
                        for i = rowstep:rowstep:xNodes
                            nd = i;
                            ndElim(nd) = nd - 1;
                        end
                        xInc = 2/obj.edgeEls(1);
                        xPrev = -1;
                        for i = 2:xNodes
                            nd = i;
                            if(ndElim(nd) == 0)
                                nodes(nd,1) = xPrev + xInc;
                                xPrev = xPrev + xInc;
                            end
                        end
                    elseif(obj.edgeEls(3) < obj.edgeEls(1))
                        nds2Elim = obj.edgeEls(1) - obj.edgeEls(3);
                        rowstep = ceil((xNodes-nds2Elim)/(nds2Elim+1)) + 1;
                        for i = rowstep:rowstep:xNodes
                            nd = (yNodes-1)*xNodes + i;
                            ndElim(nd) = nd - 1;
                        end
                        xInc = 2/obj.edgeEls(3);
                        xPrev = -1;
                        for i = 2:xNodes
                            nd = (yNodes-1)*xNodes + i;
                            if(ndElim(nd) == 0)
                                nodes(nd,1) = xPrev + xInc;
                                xPrev = xPrev + xInc;
                            end
                        end
                    end
                    if(obj.edgeEls(2) < obj.edgeEls(4))
                        nds2Elim = obj.edgeEls(4) - obj.edgeEls(2);
                        rowstep = ceil((yNodes-nds2Elim)/(nds2Elim+1)) + 1;
                        for j = rowstep:rowstep:yNodes
                            nd = j*xNodes;
                            ndElim(nd) = nd - xNodes;
                        end
                        yInc = 2/obj.edgeEls(2);
                        yPrev = -1;
                        for j = 2:yNodes
                            nd = j*xNodes;
                            if(ndElim(nd) == 0)
                                nodes(nd,2) = yPrev + yInc;
                                yPrev = yPrev + yInc;
                            end
                        end
                    elseif(obj.edgeEls(4) < obj.edgeEls(2))
                        nds2Elim = obj.edgeEls(2) - obj.edgeEls(4);
                        rowstep = ceil((yNodes-nds2Elim)/(nds2Elim+1)) + 1;
                        for j = rowstep:rowstep:yNodes
                            nd = (j-1)*xNodes + 1;
                            ndElim(nd) = nd - xNodes;
                        end
                        yInc = 2/obj.edgeEls(4);
                        yPrev = -1;
                        for j = 2:yNodes
                            nd = (j-1)*xNodes + 1;
                            if(ndElim(nd) == 0)
                                nodes(nd,2) = yPrev + yInc;
                                yPrev = yPrev + yInc;
                            end
                        end
                    end
                    if(sum(ndElim,'all') > 0.1)
                        newNds = [];
                        newLabel = zeros(length(nodes),1);
                        lab = 0;
                        for i = 1:length(nodes)
                            if(ndElim(i) == 0)
                                lab = lab + 1;
                                newLabel(i) = lab;
                                newNds = [newNds;nodes(i,:)];
                            end
                        end
                        nodes = newNds;
                        for i = 1:length(elements)
                            newEl = [];
                            for j = 1:4
                                orig = elements(i,j);
                                if(orig ~= 0)
                                    if(ndElim(orig) ~= 0)
                                        nlab = newLabel(ndElim(orig));
                                    else
                                        nlab = newLabel(orig);
                                    end
                                    if(nlab ~= 0)
                                        if(~any(newEl == nlab))
                                            newEl = [newEl,nlab];
                                        end
                                    end
                                end
                            end
                            ln = length(newEl);
                            if(ln < 4)
                                newEl = [newEl,zeros(1,4-ln)];
                            end
                            elements(i,:) = newEl;
                        end
                    end
                    newNodes = [];
                    for i = 1:length(nodes)
                        XYZ = obj.XYZCoord(nodes(i,:));
                        newNodes = [newNodes;XYZ];
                    end
                    nodes = newNodes;
                else
                    error('Only quadrilateral or cylinder shell regions can use the structured meshing option')
                end
            else
                [boundaryNodes,edges] = obj.initialBoundary();
%                 plot2DMesh(boundaryNodes,[]);
%                 keyboard
                mesh = NuMesh2D(boundaryNodes,edges(:,1:2));
                [etaNodes,elements,mesh] = mesh.createPlanarMesh(elType,1);
                nodes = [];
                for i = 1:length(etaNodes)
                    XYZ = obj.XYZCoord(etaNodes(i,:));
                    nodes = [nodes;XYZ];
                end
            end
        end
        
        function [nodes,edges] = initialBoundary(obj)
            if(contains(obj.type,'quad'))
                nodes = [];
                edges = [];
                %% initialize nodes on side 1
                delE = 2/obj.edgeEls(1);
                for i = 1:obj.edgeEls(1)
                    e1 = delE*i-1;
                    pt = [e1,-1];
                    nodes = [nodes;pt];
                end
                %% initialize nodes on side 2
                delE = 2/obj.edgeEls(2);
                for i = 1:obj.edgeEls(2)
                    e2 = delE*i-1;
                    pt = [1,e2];
                    nodes = [nodes;pt];
                end
                %% initialize nodes on side 3
                delE = 2/obj.edgeEls(3);
                for i = 1:obj.edgeEls(3)
                    e1 = 1-delE*i;
                    pt = [e1,1];
                    nodes = [nodes;pt];
                end
                %% initialize nodes on side 4
                delE = 2/obj.edgeEls(4);
                for i = 1:obj.edgeEls(4)
                    e2 = 1 - delE*i;
                    pt = [-1,e2];
                    nodes = [nodes;pt];
                end
                %% create edges
                numNds = length(nodes);
                edges = [[1:numNds]',[2:numNds,1]',-ones(numNds,1),zeros(numNds,3)];
            elseif(contains(obj.type,'tri'))
                nodes = [];
                %% initialize nodes on side 1
                delE = 1/obj.edgeEls(1);
                for i = 1:obj.edgeEls(1)
                    e1 = delE*i;
                    pt = [e1,0];
                    nodes = [nodes;pt];
                end
                %% initialize nodes on side 2
                delE = 1/obj.edgeEls(2);
                for i = 1:obj.edgeEls(2)
                    e1 = 1 - delE*i;
                    e2 = delE*i;
                    pt = [e1,e2];
                    nodes = [nodes;pt];
                end
                %% initialize nodes on side 3
                delE = 1/obj.edgeEls(3);
                for i = 1:obj.edgeEls(3)
                    e2 = 1 - delE*i;
                    pt = [0,e2];
                    nodes = [nodes;pt];
                end
                %% create edges
                numNds = length(nodes);
                edges = [[1:numNds]',[2:numNds,1]',-ones(numNds,1),zeros(numNds,3)];
            elseif(contains(obj.type,'sphere'))
                nodes = [];
                delTheta = 2*pi/obj.edgeEls(1);
                for i = 1:obj.edgeEls(1)
                    theta = delTheta*i;
                    x = 0.5*pi*cos(theta);
                    y = 0.5*pi*sin(theta);
                    nodes = [nodes;[x,y]];
                end
                numNds = length(nodes);
                edges = [[1:numNds]',[2:numNds,1]',-ones(numNds,1),zeros(numNds,3)];
            end
        end
        
        function [XYZ] = XYZCoord(obj,eta)
            switch obj.type
                case 'quad4'
                    Nvec = zeros(1,4);
                    Nvec(1) = 0.25*(eta(1)-1)*(eta(2)-1);
                    Nvec(2) = -0.25*(eta(1)+1)*(eta(2)-1);
                    Nvec(3) = 0.25*(eta(1)+1)*(eta(2)+1);
                    Nvec(4) = -0.25*(eta(1)-1)*(eta(2)+1);
                    XYZ = Nvec*obj.keyPts;
                case 'quad9'
                    r1 = -1;
                    r2 = 0;
                    r3 = 1;
                    Nvec = zeros(1,9);
                    Nvec(1) = 0.25*(eta(1)-r2)*(eta(1)-r3)*(eta(2)-r2)*(eta(2)-r3);
                    Nvec(2) = 0.25*(eta(1)-r1)*(eta(1)-r2)*(eta(2)-r2)*(eta(2)-r3);
                    Nvec(3) = 0.25*(eta(1)-r1)*(eta(1)-r2)*(eta(2)-r1)*(eta(2)-r2);
                    Nvec(4) = 0.25*(eta(1)-r2)*(eta(1)-r3)*(eta(2)-r1)*(eta(2)-r2);
                    Nvec(5) = -0.5*(eta(1)-r1)*(eta(1)-r3)*(eta(2)-r2)*(eta(2)-r3);
                    Nvec(6) = -0.5*(eta(1)-r1)*(eta(1)-r2)*(eta(2)-r1)*(eta(2)-r3);
                    Nvec(7) = -0.5*(eta(1)-r1)*(eta(1)-r3)*(eta(2)-r1)*(eta(2)-r2);
                    Nvec(8) = -0.5*(eta(1)-r2)*(eta(1)-r3)*(eta(2)-r1)*(eta(2)-r3);
                    Nvec(9) = (eta(1)-r1)*(eta(1)-r3)*(eta(2)-r1)*(eta(2)-r3);
                    XYZ = Nvec*obj.keyPts;
                case 'quad16'
                    r1 = -1;
                    r2 = -0.333333333333333;
                    r3 = 0.333333333333333;
                    r4 = 1;
                    coef = [0.31640625,-0.31640625,0.31640625,-0.31640625,...
                      -0.94921875,0.94921875,0.94921875,-0.94921875,...
                      -0.94921875,0.94921875,0.94921875,-0.94921875,...
                       2.84765625,-2.84765625,2.84765625,-2.84765625];
                    Nvec = zeros(1,16);
                    Nvec(1) = (eta(1)-r2)*(eta(1)-r3)*(eta(1)-r4)*(eta(2)-r2)*(eta(2)-r3)*(eta(2)-r4);
                    Nvec(2) = (eta(1)-r1)*(eta(1)-r2)*(eta(1)-r3)*(eta(2)-r2)*(eta(2)-r3)*(eta(2)-r4);
                    Nvec(3) = (eta(1)-r1)*(eta(1)-r2)*(eta(1)-r3)*(eta(2)-r1)*(eta(2)-r2)*(eta(2)-r3);
                    Nvec(4) = (eta(1)-r2)*(eta(1)-r3)*(eta(1)-r4)*(eta(2)-r1)*(eta(2)-r2)*(eta(2)-r3);
                    Nvec(5) = (eta(1)-r1)*(eta(1)-r3)*(eta(1)-r4)*(eta(2)-r2)*(eta(2)-r3)*(eta(2)-r4);
                    Nvec(6) = (eta(1)-r1)*(eta(1)-r2)*(eta(1)-r4)*(eta(2)-r2)*(eta(2)-r3)*(eta(2)-r4);
                    Nvec(7) = (eta(1)-r1)*(eta(1)-r2)*(eta(1)-r3)*(eta(2)-r1)*(eta(2)-r3)*(eta(2)-r4);
                    Nvec(8) = (eta(1)-r1)*(eta(1)-r2)*(eta(1)-r3)*(eta(2)-r1)*(eta(2)-r2)*(eta(2)-r4);
                    Nvec(9) = (eta(1)-r1)*(eta(1)-r2)*(eta(1)-r4)*(eta(2)-r1)*(eta(2)-r2)*(eta(2)-r3);
                    Nvec(10) = (eta(1)-r1)*(eta(1)-r3)*(eta(1)-r4)*(eta(2)-r1)*(eta(2)-r2)*(eta(2)-r3);
                    Nvec(11) = (eta(1)-r2)*(eta(1)-r3)*(eta(1)-r4)*(eta(2)-r1)*(eta(2)-r2)*(eta(2)-r4);
                    Nvec(12) = (eta(1)-r2)*(eta(1)-r3)*(eta(1)-r4)*(eta(2)-r1)*(eta(2)-r3)*(eta(2)-r4);
                    Nvec(13) = (eta(1)-r1)*(eta(1)-r3)*(eta(1)-r4)*(eta(2)-r1)*(eta(2)-r3)*(eta(2)-r4);
                    Nvec(14) = (eta(1)-r1)*(eta(1)-r2)*(eta(1)-r4)*(eta(2)-r1)*(eta(2)-r3)*(eta(2)-r4);
                    Nvec(15) = (eta(1)-r1)*(eta(1)-r2)*(eta(1)-r4)*(eta(2)-r1)*(eta(2)-r2)*(eta(2)-r4);
                    Nvec(16) = (eta(1)-r1)*(eta(1)-r3)*(eta(1)-r4)*(eta(2)-r1)*(eta(2)-r2)*(eta(2)-r4);
                    Nvec = coef.*Nvec;
                    XYZ = Nvec*obj.keyPts;
                case 'tri3'
                    Nvec = zeros(1,3);
                    Nvec(1) = 1 - eta(1) - eta(2);
                    Nvec(2) = eta(1);
                    Nvec(3) = eta(2);
                    XYZ = Nvec*obj.keyPts;
                case 'tri6'
                    Nvec = zeros(1,6);
                    Nvec(1) = 2*(eta(1) + eta(2) - 1)*(eta(1) + eta(2) - 0.5);
                    Nvec(2) = 2*eta(1)*(eta(1)-0.5);
                    Nvec(3) = 2*eta(2)*(eta(2)-0.5);
                    Nvec(4) = -4*eta(1)*(eta(1) + eta(2) - 1);
                    Nvec(5) = 4*eta(1)*eta(2);
                    Nvec(6) = -4*eta(2)*(eta(1) + eta(2) - 1);
                    XYZ = Nvec*obj.keyPts;
                case 'tri10'
                    r2 = 1/3;
                    r3 = 2/3;
                    coef = [-4.5,4.5,4.5,13.5,-13.5,13.5,13.5,-13.5,13.5,-27];
                    Nvec = zeros(1,10);
                    Nvec(1) = (eta(1)+eta(2)-r2)*(eta(1)+eta(2)-r3)*(eta(1)+eta(2)-1);
                    Nvec(2) = eta(1)*(eta(1)-r2)*(eta(1)-r3);
                    Nvec(3) = eta(2)*(eta(2)-r2)*(eta(2)-r3);
                    Nvec(4) = eta(1)*(eta(1)+eta(2)-r3)*(eta(1)+eta(2)-1);
                    Nvec(5) = eta(1)*(eta(1)-r2)*(eta(1)+eta(2)-1);
                    Nvec(6) = eta(1)*eta(2)*(eta(1)-r2);
                    Nvec(7) = eta(1)*eta(2)*(eta(2)-r2);
                    Nvec(8) = eta(2)*(eta(2)-r2)*(eta(1)+eta(2)-1);
                    Nvec(9) = eta(2)*(eta(1)+eta(2)-r3)*(eta(1)+eta(2)-1);
                    Nvec(10) = eta(1)*eta(2)*(eta(1)+eta(2)-1);
                    Nvec = coef.*Nvec;
                    XYZ = Nvec*obj.keyPts;
                case 'sphere'
                    vec = obj.keyPts(2,:) - obj.keyPts(1,:);
                    outerRad = sqrt(vec*vec');
                    phiComp = sqrt(eta*eta');
                    if(eta(1) > 0)
                        theta = atan(eta(2)/eta(1));
                    elseif(eta(1) < 0)
                        theta = pi + atan(eta(2)/eta(1));
                    else
                        theta = atan(eta(2)/1e-8);
                    end
                    phi = 0.5*pi - phiComp;
                    xloc = outerRad*cos(theta)*cos(phi);
                    yloc = outerRad*sin(theta)*cos(phi);
                    zloc = outerRad*sin(phi);
                    a1 = (1/outerRad)*vec;
                    vec2 = obj.keyPts(3,:) - obj.keyPts(1,:);
                    vec3 = [(vec(2)*vec2(3) - vec(3)*vec2(2)),...
                        (vec(3)*vec2(1) - vec(1)*vec2(3)),...
                        (vec(1)*vec2(2) - vec(2)*vec2(1))];
                    mag = sqrt(vec3*vec3');
                    a3 = (1/mag)*vec3;
                    a2 = [(a3(2)*a1(3) - a3(3)*a1(2)),...
                        (a3(3)*a1(1) - a3(1)*a1(3)),...
                        (a3(1)*a1(2) - a3(2)*a1(1))];
                    alpha = [a1;a2;a3];
                    XYZ = ([xloc,yloc,zloc]*alpha + obj.keyPts(1,:));
            end
        end

    end
end

