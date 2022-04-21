classdef NuMesh2D
    % Object representing a 2D mesh of triangular or quadrilateral elements
    %   
    
    properties
        nodes = [];
        edges = [];
        triElements = [];
        quadElements = [];
        nodeGL;
        edgeGL;
        triElGL;
        maxEdgeLen = 1;
        avgProjLen = 1;
        maxElSize = 1;
        minElSize = 1;
    end
    
    methods
        function obj = NuMesh2D(boundaryNodes,boundaryEdges)
            obj.nodes = boundaryNodes;
            numNds = length(boundaryNodes);
            if(isempty(boundaryEdges))
                obj.edges = zeros(numNds,6);
                obj.edges(:,1) = [1:numNds]';
                obj.edges(:,2) = [2:numNds,1]';
                obj.edges(:,3) = -1;
            else
                numEdges = size(boundaryEdges,1);
                obj.edges = [boundaryEdges,-ones(numEdges,1),zeros(numEdges,3)];
            end
        end
        
        function [nodeFound,edgesAdded,obj] = adoptAnyNode(obj,currentEdge,ndPt,srchRad)
            nodeFound = 0;
            edgesAdded = 0;
            n1 = obj.edges(currentEdge,1);
            n2 = obj.edges(currentEdge,2);
            edgeMpt = 0.5*(obj.nodes(n1,:) + obj.nodes(n2,:));
            nearNds = obj.nodeGL.findInRadius(ndPt,srchRad);
            k = 1;
            while(nodeFound == 0 && k <= length(nearNds))
                j = nearNds(k);
                if(j ~= n1 && j ~= n2)
                    vec = obj.nodes(j,:) - ndPt;
                    dist = sqrt(vec*vec');
                    mpVec = obj.nodes(j,:) - edgeMpt;
                    dp = obj.edges(currentEdge,5:6)*mpVec';
                    if(dist < srchRad && dp > 0)
                        elNum = obj.findElement([n1,n2,j]);
                        if(elNum == 0)
                            [violation] = obj.checkViolations([n1,n2,j],[],1.3);
                            if(violation == 0)
                                    [numAdded,obj] = obj.completeEdges(n1,n2,j,size(obj.triElements,1)+1);
                                    edgesAdded = edgesAdded + numAdded;
                                    newEl = [n1,n2,j];
                                    obj.triElements = [obj.triElements;newEl];
                                    midPt = 0.333333*(obj.nodes(n1,:) + obj.nodes(n2,:) + obj.nodes(j,:));
                                    obj.triElGL = obj.triElGL.addEntry(size(obj.triElements,1),midPt);
                                    obj.edges(currentEdge,4) = size(obj.triElements,1);
                                    nodeFound = 1;
                            end
                        else
                            obj.edges(currentEdge,4) = elNum;
                            nodeFound = 1;
                        end
                    end
                end
                k = k + 1;
            end
        end
        
        function [nodeFound,edgesAdded,obj] = adoptConnectedNode(obj,currentEdge,ndPt,srchRad)
            nodeFound = 0;
            edgesAdded = 0;
            n1 = obj.edges(currentEdge,1);
            n2 = obj.edges(currentEdge,2);
            edgeMpt = 0.5*(obj.nodes(n1,:) + obj.nodes(n2,:));
            nearEdges = obj.edgeGL.findInRadius(edgeMpt,1.1*obj.maxEdgeLen);
            k = 1;
            while(nodeFound == 0 && k <= length(nearEdges))
                j = nearEdges(k);
                if(currentEdge ~= j)
                    n1j = obj.edges(j,1);
                    n2j = obj.edges(j,2);
                    if(n1j == n1 || n1j == n2 || n2j == n1 || n2j == n2)
                        if(n1j == n1 || n1j == n2)
                            n3 = n2j;
                        else
                            n3 = n1j;
                        end
                        vec = obj.nodes(n3,1:2) - ndPt;
                        dist = sqrt(vec*vec');
                        mpVec = obj.nodes(n3,:) - edgeMpt;
                        dp = obj.edges(currentEdge,5:6)*mpVec';
                        if(dist < srchRad && dp > 0)
                            elNum = obj.findElement([n1,n2,n3]);
                            if(elNum == 0)
                                violation = obj.checkViolations([n1,n2,n3],[],1.3);
                                if(violation == 0)
                                        [numAdded,obj] = obj.completeEdges(n1,n2,n3,size(obj.triElements,1)+1);
                                        edgesAdded = edgesAdded + numAdded;
                                        newEl = [n1,n2,n3];
                                        obj.triElements = [obj.triElements;newEl];
                                        midPt = 0.333333*(obj.nodes(n1,:) + obj.nodes(n2,:) + obj.nodes(n3,:));
                                        obj.triElGL = obj.triElGL.addEntry(size(obj.triElements,1),midPt);
                                        obj.edges(currentEdge,4) = size(obj.triElements,1);
                                        nodeFound = 1;
                                end
                            else
                                obj.edges(currentEdge,4) = elNum;
                                nodeFound = 1;
                            end
                        end
                    end
                end
                k = k + 1;
            end
        end
        
        function [violation] = checkViolations(obj,ndLab,ndCrd,marginFact)
            violation = 0;
            elCrd = [];
            for i = 1:length(ndLab)
                if(ndLab(i) > 0)
                    elCrd = [elCrd;obj.nodes(ndLab(i),:)];
                else
                    j = -ndLab(i);
                    elCrd = [elCrd;ndCrd(j,:)];
                end
            end
            midPt = 0.333333*(elCrd(1,:) + elCrd(2,:) + elCrd(3,:));
            nearEls = obj.triElGL.findInRadius(midPt,2.1*obj.maxEdgeLen);
            for m = 1:length(nearEls)
                nm = obj.triElements(nearEls(m),:);
                thisElnds = [obj.nodes(nm(1),:);obj.nodes(nm(2),:);obj.nodes(nm(3),:)];
                overlap = obj.elementsOverlap([1,2,3;4,5,6],[elCrd;thisElnds],1);
                if(overlap == 1)
                    violation = 1;
                end
                for i = 1:3
                    if(ndLab(i) < 0)
                        inEl = obj.ptInEl(elCrd(i,:),[1,2,3],thisElnds,marginFact);
                        if(inEl == 1)
                            violation = 1;
                        end
                    end
                end
            end
            nearNds = obj.nodeGL.findInRadius(midPt,1.1*obj.maxEdgeLen);
            for m = 1:length(nearNds)
                nd = nearNds(m);
                if(nd ~= ndLab(1) && nd ~= ndLab(2) && nd ~= ndLab(3))
                    pt = obj.nodes(nd,:);
                    inEl = obj.ptInEl(pt,[1,2,3],elCrd,marginFact);
                    if(inEl == 1)
                        violation = 1;
                    end
                end
            end
        end

        function [numAdded,obj] = completeEdges(obj,n1,n2,n3,elNum)
            numAdded = 0;
            midPt = 0.333333*(obj.nodes(n1,:) + obj.nodes(n2,:) + obj.nodes(n3,:));
            nearEdges = obj.edgeGL.findInRadius(midPt,1.1*obj.maxEdgeLen);
            e2Found = 0;
            j = 1;
            while(e2Found == 0 && j <= length(nearEdges))
                k = nearEdges(j);
                if((obj.edges(k,1) == n1 && obj.edges(k,2) == n3) || (obj.edges(k,1) == n3 && obj.edges(k,2) == n1))
                    e2Found = k;
                end
                j = j + 1;
            end
            if(e2Found == 0)
                newEdge = [n1,n3,elNum,0,0,0];
                obj.edges = [obj.edges;newEdge];
                numAdded = numAdded + 1;
                midPt = 0.5*(obj.nodes(n1,:) + obj.nodes(n3,:));
                obj.edgeGL = obj.edgeGL.addEntry(length(obj.edges),midPt);
                vec = obj.nodes(n1,:) - obj.nodes(n3,:);
                mag = sqrt(vec*vec');
                if(mag > obj.maxEdgeLen)
                    obj.maxEdgeLen = mag;
                end
            else
                obj.edges(e2Found,4) = elNum;
            end
            e3Found = 0;
            j = 1;
            while(e3Found == 0 && j <= length(nearEdges))
                k = nearEdges(j);
                if((obj.edges(k,1) == n2 && obj.edges(k,2) == n3) || (obj.edges(k,1) == n3 && obj.edges(k,2) == n2))
                    e3Found = k;
                end
                j = j + 1;
            end
            if(e3Found == 0)
                newEdge = [n2,n3,elNum,0,0,0];
                obj.edges = [obj.edges;newEdge];
                numAdded = numAdded + 1;
                midPt = 0.5*(obj.nodes(n2,:) + obj.nodes(n3,:));
                obj.edgeGL = obj.edgeGL.addEntry(length(obj.edges),midPt);
                vec = obj.nodes(n2,:) - obj.nodes(n3,:);
                mag = sqrt(vec*vec');
                if(mag > obj.maxEdgeLen)
                    obj.maxEdgeLen = mag;
                end
            else
                obj.edges(e3Found,4) = elNum;
            end
        end
        
        function [nodeCreated,edgesAdded,obj] = createNode(obj,currentEdge,ndPt)
            nodeCreated = 0;
            edgesAdded = 0;
            n1 = obj.edges(currentEdge,1);
            n2 = obj.edges(currentEdge,2);
            violation = obj.checkViolations([n1,n2,-1],ndPt,1.3);
            if(violation == 0)
                newNd = ndPt;
                obj.nodes = [obj.nodes;newNd];
                obj.nodeGL = obj.nodeGL.addEntry(length(obj.nodes),newNd);
                n3 = length(obj.nodes);
                newEl = [n1,n2,n3];
                obj.triElements = [obj.triElements;newEl];
                midPt = 0.333333*(obj.nodes(n1,:) + obj.nodes(n2,:) + obj.nodes(n3,:));
                obj.triElGL = obj.triElGL.addEntry(size(obj.triElements,1),midPt);
                numEls = size(obj.triElements,1);
                [numAdded,obj] = obj.completeEdges(n1,n2,n3,numEls);
                edgesAdded = edgesAdded + numAdded;
                nodeCreated = 1;
            end
        end
        
        function [nodes,elements,obj] = createPlanarMesh(obj,elType,equalizeSpacing)
            if(contains(elType,'quad'))
                obj = obj.skewNodes();
            end
            boundaryNodes = obj.nodes;
            if(equalizeSpacing)
                maxIt = ceil(0.5*length(obj.nodes));
                obj = obj.uniformBoundarySpacing(maxIt);
            end
%             plot2DMesh(obj.nodes,[]);
%             keyboard
            obj = obj.createUnstructTriMesh();
%             plot2DMesh(obj.nodes,obj.triElements);
%             keyboard
            obj = obj.distributeNodes(boundaryNodes);
%             plot2DMesh(obj.nodes,obj.triElements);
%             keyboard
            if(contains(elType,'quad'))
                obj = obj.unSkewNodes();
%                 plot2DMesh(obj.nodes,obj.triElements);
%                 keyboard
            end
            obj = obj.mergeTriEls(elType);
            nodes = obj.nodes;
            if(contains(elType,'quad'))
                elements = obj.quadElements;
            else
                elements = obj.triElements;
            end
            plot2DMesh(nodes,elements);
            keyboard;
        end
        
        function [nodes,elements,obj] = createSweptMesh(obj,sweepMethod,direction,sweepDistance,sweepElements,followNormal,destNodes)
            numNds = size(obj.nodes,1);
            numEdges = size(obj.edges,1);
            obj.quadElements = [];
            if(followNormal == 1)
                methString = 'in_direction, to_point, from_point'; 
                if(contains(methString,sweepMethod))
                    
                end
            else
                switch sweepMethod
                    case 'in_direction'
                        dirMag = sqrt(direction*direction');
                        unitDir = (1/dirMag)*direction;
                        ndDir = (sweepDistance/sweepElements)*ones(numNds,2)*[unitDir(1),0;0,unitDir(2)];
                end
            end
            if(contains(sweepMethod,'revolve'))
            else
                for i = 1:sweepElements
                    ndRow = obj.nodes(1:numNds,:) + i*ndDir;
                    obj.nodes = [obj.nodes;ndRow];
                    for j = 1:numEdges
                        n1 = (i-1)*numNds + obj.edges(j,1);
                        n2 = (i-1)*numNds + obj.edges(j,2);
                        n3 = i*numNds + obj.edges(j,2);
                        n4 = i*numNds + obj.edges(j,1);
                        obj.quadElements = [obj.quadElements;[n1,n2,n3,n4]];
                    end
                end
            end
            nodes = obj.nodes;
            elements = obj.quadElements;
        end
        
        function obj = createUnstructTriMesh(obj)
            obj = obj.prepareForMesh();
            numEdges = length(obj.edges);
            for i = 1:numEdges
                n1i = obj.edges(i,1);
                n2i = obj.edges(i,2);
                midPt = 0.5*(obj.nodes(n1i,:) + obj.nodes(n2i,:));
                nearEdges = obj.edgeGL.findInRadius(midPt,1.1*obj.maxEdgeLen);
                for m = 1:length(nearEdges)
                    j = nearEdges(m);
                    if(j ~= i)
                        n1j = obj.edges(j,1);
                        n2j = obj.edges(j,2);
                        allNds = [n1i,n2i,n1j,n2j];
                        allNds = sort(allNds);
                        n1 = 0;
                        for k = 1:3
                            if(allNds(k+1) == allNds(k))
                                n1 = allNds(k);
                                allNds(k) = 0;
                                allNds(k+1) = 0;
                            end
                        end
                        allNds = sort(allNds);
                        n2 = allNds(3);
                        n3 = allNds(4);
                        if(n1 ~= 0)
                            v1 = obj.nodes(n2,:) - obj.nodes(n1,:);
                            mag1 = sqrt(v1*v1');
                            v2 = obj.nodes(n3,:) - obj.nodes(n1,:);
                            mag2 = sqrt(v2*v2');
                            dp = (v1*v2')/(mag1*mag2);
                            if(dp > 0.4)
                                v3 = v1 + v2;
                                dp = obj.edges(i,5:6)*v3';
                                if(dp > 0)
                                    newEl = [n1,n2,n3];
                                    eNum = obj.findElement(newEl);
                                    if(eNum == 0)
                                        obj.triElements = [obj.triElements;newEl];
                                        midPt = 0.3333333*(obj.nodes(n1,:) + obj.nodes(n2,:) + obj.nodes(n3,:));
                                        obj.triElGL = obj.triElGL.addEntry(size(obj.triElements,1),midPt);
                                        eNum = size(obj.triElements,1);
                                        obj.edges(i,4) = eNum;
                                        obj.edges(j,4) = eNum;
                                        newEdge = [n2,n3,eNum,0,0,0];
                                        obj.edges = [obj.edges;newEdge];
                                        midPt = 0.5*(obj.nodes(n2,:) + obj.nodes(n3,:));
                                        obj.edgeGL = obj.edgeGL.addEntry(length(obj.edges),midPt);
                                        vec = obj.nodes(n2,:) - obj.nodes(n3,:);
                                        mag = sqrt(vec*vec');
                                        if(mag > obj.maxEdgeLen)
                                            obj.maxEdgeLen = mag;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            edgesAdded = length(obj.edges);
            maxIt = edgesAdded;
            it = 0;
            while(edgesAdded > 0 && it < maxIt)
                edgesAdded = 0;
                obj = obj.sortIncompleteEdges();
                numEdges = length(obj.edges);
%                 plot2DMesh(obj.nodes,obj.triElements);
%                 keyboard
                for i = 1:numEdges
                    if(obj.edges(i,3) > 0)
                        n1 = obj.edges(i,1);
                        n2 = obj.edges(i,2);
                        midpt = 0.5*(obj.nodes(n2,:) + obj.nodes(n1,:));
                        eVec = obj.nodes(n2,:) - obj.nodes(n1,:);
                        eNorm = [-eVec(2),eVec(1)];
                        el = obj.edges(i,3);
                        for j = 1:3
                            if(obj.triElements(el,j) ~= n1 && obj.triElements(el,j) ~= n2)
                                n3 = obj.triElements(el,j);
                                vin = midpt - obj.nodes(n3,:);
                                dp = eNorm*vin';
                                if(dp < 0)
                                    eNorm = -eNorm;
                                end
                                mag = sqrt(eNorm*eNorm');
                                obj.edges(i,5:6) = (1/mag)*eNorm;
                            end
                        end
                    end
                end
                for i = 1:numEdges
                    if(obj.edges(i,4) == 0)
                        n1 = obj.edges(i,1);
                        n2 = obj.edges(i,2);
                        midpt = 0.5*(obj.nodes(n1,1:2) + obj.nodes(n2,1:2));
                        vec = obj.nodes(n1,1:2) - obj.nodes(n2,1:2);
                        eLen = sqrt(vec*vec');
                        unitProj = obj.edges(i,5:6);
                        projLen = 0.8*obj.avgProjLen + 0.2*0.866025403784439*eLen;
                        srchRad = 0.5*sqrt(projLen^2 + eLen^2);
                        ndPt = midpt + 0.5*projLen*unitProj;
                        [nodeFound,edAdd,obj] = obj.adoptConnectedNode(i,ndPt,srchRad);
                        edgesAdded = edgesAdded + edAdd;
                        if(nodeFound == 0)
                            ndPt = midpt + projLen*unitProj;
                            [nodeFound,edAdd,obj] = obj.adoptConnectedNode(i,ndPt,srchRad);
                            edgesAdded = edgesAdded + edAdd;
                        end
                        if(nodeFound == 0)
                            ndPt = midpt + 0.5*projLen*unitProj;
                            [nodeFound,edAdd,obj] = obj.adoptAnyNode(i,ndPt,srchRad);
                            edgesAdded = edgesAdded + edAdd;
                        end
                        if(nodeFound == 0)
                            ndPt = midpt + projLen*unitProj;
                            [nodeFound,edAdd,obj] = obj.adoptAnyNode(i,ndPt,srchRad);
                            edgesAdded = edgesAdded + edAdd;
                        end
                        if(nodeFound == 0)
                            ndPt = midpt + projLen*unitProj;
                            [nodeCreated,edAdd,obj] = obj.createNode(i,ndPt);
                            edgesAdded = edgesAdded + edAdd;
                            if(nodeCreated == 0)
                                ndPt = midpt + 0.5*projLen*unitProj;
                                [nodeCreated,edAdd,obj] = obj.createNode(i,ndPt);
                                edgesAdded = edgesAdded + edAdd;
                            end
                        end
                    end
                end
%                 tic
%                 for i = 1:length(obj.edges)
%                     if(obj.edges(i,4) == 0)
%                         n1 = obj.edges(i,1);
%                         n2 = obj.edges(i,2);
%                         midPt = 0.5*(obj.nodes(n1,:) + obj.nodes(n2,:));
%                         nearEdges = obj.edgeGL.findInRadius(midPt,1.1*obj.maxEdgeLen);
%                         for p = 1:length(nearEdges)
%                             j = nearEdges(p);
%                             if(j ~= i && obj.edges(j,4) == 0 && obj.edges(i,4) == 0)
%                                 for q = 1:length(nearEdges)
%                                     k = nearEdges(q);
%                                     if(k ~= i && k ~= j && obj.edges(i,4) == 0 && obj.edges(j,4) == 0 && obj.edges(k,4) == 0)
%                                         allNds = [obj.edges(i,1:2),obj.edges(j,1:2),obj.edges(k,1:2)];
%                                         allNds = sort(allNds);
%                                         numDup = 0;
%                                         for m = 1:5
%                                             if(allNds(m+1) == allNds(m))
%                                                 numDup = numDup + 1;
%                                             end
%                                         end
%                                         if(numDup == 3)
%                                             newEl = [allNds(1),allNds(3),allNds(5)];
%                                             elNum = obj.findElement(newEl);
%                                             if(elNum == 0)
%                                                 violation = obj.checkViolations(newEl,[],1.3);
%                                                 if(violation == 0)
%                                                     obj.triElements = [obj.triElements;newEl];
%                                                     elNum = size(obj.triElements,1);
%                                                     midPt = 0.333333*(obj.nodes(newEl(1),:) + obj.nodes(newEl(2),:) + obj.nodes(newEl(3),:));
%                                                     obj.triElGL = obj.triElGL.addEntry(elNum,midPt);
%                                                     obj.edges(i,4) = elNum;
%                                                     obj.edges(j,4) = elNum;
%                                                     obj.edges(k,4) = elNum;
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%                 toc
                it = it + 1;
            end
        end
        
        function obj = distributeNodes(obj,boundaryNodes)
            numNds = length(obj.nodes);
            dim = 2*numNds;
            Dmat = zeros(dim,1);
            numBound = size(boundaryNodes,1);
            Dmat(1:2*numBound) = 100000;
            Pmat = 10*ones(dim,1) + Dmat;
            Pinv = 1./Pmat;
            rhs = zeros(dim,1);
            for i = 1:numBound
                j = 2*i;
                rhs(j-1:j) = 100000*boundaryNodes(i,:)';
            end
            
            numEls = size(obj.triElements,1);
            elWt = zeros(numEls,1);
            for i = 1:numEls
                ni = obj.triElements(i,:);
                v1 = obj.nodes(ni(2),:) - obj.nodes(ni(1),:);
                v2 = obj.nodes(ni(3),:) - obj.nodes(ni(1),:);
                cp = v1(1)*v2(2) - v1(2)*v2(1);
                elWt(i) = abs(cp);
            end
            avgWt = mean(elWt,'all');
            elWt = (1/avgWt)*elWt;

            eMat = [2     0    -1     0    -1     0;...
                    0     2     0    -1     0    -1;...
                   -1     0     2     0    -1     0;...
                    0    -1     0     2     0    -1;...
                   -1     0    -1     0     2     0;...
                    0    -1     0    -1     0     2];

            xVec = zeros(dim,1);
            gVec = -rhs;
            wVec = Pinv.*gVec;
            hVec = -wVec;
            res = gVec'*wVec;
            i = 0;
            while(res > 1e-12 && i < dim)
                zVec = zeros(dim,1);
                for j=1:size(obj.triElements,1)
                    for k = 1:3
                        nd = obj.triElements(j,k);
                        for m = 1:2
                            ei = 2*(k-1) + m;
                            gi = 2*(nd-1) + m;
                            for n = 1:3
                                ndj = obj.triElements(j,n);
                                for p = 1:2
                                    ej = 2*(n-1) + p;
                                    gj = 2*(ndj-1) + p;
                                    zVec(gi) = zVec(gi) + elWt(j)*eMat(ei,ej)*hVec(gj);
                                end
                            end
                        end
                    end
                end
                zVec = zVec + Dmat.*hVec;
                alpha = res/(hVec'*zVec);
                xVec = xVec + alpha*hVec;
                gVec = gVec + alpha*zVec;
                wVec = Pinv.*gVec;
                rNext = gVec'*wVec;
                beta = rNext/res;
                res = rNext;
                hVec = -wVec + beta*hVec;
                i = i + 1;
            end
            for i = 1:numNds
                j = i*2;
                obj.nodes(i,:) = xVec(j-1:j)'; 
            end
        end
        
        function [overlap] = elementsOverlap(obj,els,ndCrd,marginFact)
            Xel1 = [ndCrd(els(1,1),:);ndCrd(els(1,2),:);ndCrd(els(1,3),:)];
            Xel2 = [ndCrd(els(2,1),:);ndCrd(els(2,2),:);ndCrd(els(2,3),:)];
            cent = 0.3333333*(Xel1(1,:) + Xel1(2,:) + Xel1(3,:));
            centMat = ones(3,2)*[cent(1),0;0,cent(2)];
            Xel1 = marginFact*(Xel1 - centMat) + centMat;
            elEdges = [1,2;2,3;3,1];
            overlap = 0;
            for i = 1:3
                n11 = elEdges(i,1);
                n21 = elEdges(i,2);
                x01 = Xel1(n11,:)';
                v1 = (Xel1(n21,:) - Xel1(n11,:))';
                for j = 1:3
                    n12 = elEdges(j,1);
                    n22 = elEdges(j,2);
                    x02 = Xel2(n12,:)';
                    v2 = (Xel2(n22,:) - Xel2(n12,:))';
                    vMat = [v1,-v2];
                    mag = max(abs(vMat),[],'all');
                    [Q,R] = qr(vMat,0);
                    if(sqrt(abs(R(1,1)*R(2,2))) > 1e-6*mag)
                        soln = R\Q'*(x02 - x01);
                        if(soln(1) > 0.000001 && soln(1) < 0.999999 && soln(2) > 0.000001 && soln(2) < 0.999999)
                            overlap = 1;
                        end
                    end
                end
            end
        end

        function [elNum] = findElement(obj,nds)
            nearEls = obj.triElGL.findInRadius(obj.nodes(nds(1),:),1.1*obj.maxEdgeLen);
            sortedNds = sort(nds);
            elNum = 0;
            for j = 1:length(nearEls)
                i = nearEls(j);
                sortedEl = sort(obj.triElements(i,:));
                if(sortedNds == sortedEl)
                    elNum = i;
                end
            end
        end
        
        function [edgeDir] = getBoundaryData(obj)
            numEdges = length(obj.edges);
            avgSpacing = 0;
            edgeDir = [];
            for i = 1:numEdges
                n1 = obj.edges(i,1);
                n2 = obj.edges(i,2);
                vec = obj.nodes(n2,:) - obj.nodes(n1,:);
                mag = sqrt(vec*vec');
                avgSpacing = avgSpacing + mag;
                unitE = (1/mag)*vec;
                edgeDir = [edgeDir;unitE];
            end

            avgSpacing = avgSpacing/numEdges;
            edgeDir = avgSpacing*edgeDir;
        end
        
        function obj = getBoundaryEdgeNormals(obj,Xmin,Xmax,Ymin,Ymax,spacing)
            for i = 1:size(obj.edges,1)
                n1 = obj.edges(i,1);
                n2 = obj.edges(i,2);
                v2 = obj.nodes(n2,:) - obj.nodes(n1,:);
                eNorm = [-v2(2),v2(1)];
                mag = sqrt(eNorm*eNorm');
                obj.edges(i,5:6) = (1/mag)*eNorm;
            end
            for x = Xmin:spacing:Xmax
                v1 = [0;1];
                x01 = [x;Ymin];
                intersects = [];
                nearEdges = obj.edgeGL.findInXYMargin(x01,0.6*obj.maxEdgeLen,-1);
                for j = 1:length(nearEdges)
                    i = nearEdges(j);
                    n1 = obj.edges(i,1);
                    n2 = obj.edges(i,2);
                    v2 = obj.nodes(n2,:)' - obj.nodes(n1,:)';
                    x02 = obj.nodes(n1,:)';
                    vMat = [v1,-v2];
                    mag = max(abs(vMat),[],'all');
                    [Q,R] = qr(vMat,0);
                    if(sqrt(abs(R(1,1)*R(2,2))) > 1e-6*mag)
                        soln = R\Q'*(x02-x01);
                        if(soln(2) > 0 && soln(2) < 1)
                            intersects = [intersects;[i,soln(1)]];
                        end
                    end
                end
                numInt = size(intersects,1);
                if(abs(round(0.5*numInt) - 0.5*numInt) < 0.001)
                    for i = 1:numInt
                        for j = 1:numInt-1
                            if(intersects(j+1,2) < intersects(j,2))
                                swap = intersects(j+1,:);
                                intersects(j+1,:) = intersects(j,:);
                                intersects(j,:) = swap;
                            end
                        end
                    end
                    for i = 1:2:numInt
                        ed = intersects(i,1);
                        if(obj.edges(ed,6) < 0)
                            obj.edges(ed,5:6) = -obj.edges(ed,5:6);
                        end
                    end
                    for i = 2:2:numInt
                        ed = intersects(i,1);
                        if(obj.edges(ed,6) > 0)
                            obj.edges(ed,5:6) = -obj.edges(ed,5:6);
                        end
                    end
                end
            end
            for y = Ymin:spacing:Ymax
                v1 = [1;0];
                x01 = [Xmin;y];
                intersects = [];
                nearEdges = obj.edgeGL.findInXYMargin(x01,-1,0.6*obj.maxEdgeLen);
                for j = 1:length(nearEdges)
                    i = nearEdges(j);
                    n1 = obj.edges(i,1);
                    n2 = obj.edges(i,2);
                    v2 = obj.nodes(n2,:)' - obj.nodes(n1,:)';
                    x02 = obj.nodes(n1,:)';
                    vMat = [v1,-v2];
                    mag = max(abs(vMat),[],'all');
                    [Q,R] = qr(vMat,0);
                    if(sqrt(abs(R(1,1)*R(2,2))) > 1e-6*mag)
                        soln = R\Q'*(x02-x01);
                        if(soln(2) > 0 && soln(2) < 1)
                            intersects = [intersects;[i,soln(1)]];
                        end
                    end
                end
                numInt = size(intersects,1);
                if(abs(round(0.5*numInt) - 0.5*numInt) < 0.001)
                    for i = 1:numInt
                        for j = 1:numInt-1
                            if(intersects(j+1,2) < intersects(j,2))
                                swap = intersects(j+1,:);
                                intersects(j+1,:) = intersects(j,:);
                                intersects(j,:) = swap;
                            end
                        end
                    end
                    for i = 1:2:numInt
                        ed = intersects(i,1);
                        if(obj.edges(ed,5) < 0)
                            obj.edges(ed,5:6) = -obj.edges(ed,5:6);
                        end
                    end
                    for i = 2:2:numInt
                        ed = intersects(i,1);
                        if(obj.edges(ed,5) > 0)
                            obj.edges(ed,5:6) = -obj.edges(ed,5:6);
                        end
                    end
                end
            end
        end
        
        function [newElMerged,obj] = mergeElsBetweenAngles(obj,elMerged,el2El,minAngle,maxAngle)
            for i = 1:size(obj.triElements,1)
                if(elMerged(i) == 0)
                    for p = 1:3
                        j = el2El(i,p);
                        if(j > 0)
                            if(elMerged(i) == 0 && elMerged(j) == 0 && i ~= j)
                                allNds = [obj.triElements(i,:),obj.triElements(j,:)];
                                allNds = sort(allNds);
                                commonNds = [];
                                for k = 1:5
                                    if(allNds(k) == allNds(k+1))
                                        commonNds = [commonNds,allNds(k)];
                                        allNds(k) = 0;
                                        allNds(k+1) = 0;
                                    end
                                end
                                if(length(commonNds) == 2)
                                    allNds = sort(allNds);
                                    unCommon = allNds(5:6);
                                    quadNds = [unCommon(1),commonNds(1),unCommon(2),commonNds(2)];
                                    next = [quadNds(2:4),quadNds(1)];
                                    last = [quadNds(4),quadNds(1:3)];
                                    mergeOk = 1;
                                    for k = 1:4
                                        v1 = obj.nodes(next(k),:) - obj.nodes(quadNds(k),:);
                                        mag1 = sqrt(v1*v1');
                                        v2 = obj.nodes(last(k),:) - obj.nodes(quadNds(k),:);
                                        mag2 = sqrt(v2*v2');
                                        dp = v1*v2';
                                        theta = acos(dp/(mag1*mag2));
                                        if(theta < minAngle || theta > maxAngle)
                                            mergeOk = 0;
                                        end
                                    end
                                    if(mergeOk == 1)
                                        obj.quadElements = [obj.quadElements;quadNds];
                                        elMerged(i) = 1;
                                        elMerged(j) = 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            newElMerged = elMerged;
        end
        
        function [newElMerged,newNodeElim,obj] = mergeNestedQuadEls(obj,elMerged,el2El,nodeElim)
            numEls = size(obj.triElements,1);
            for i = 1:numEls
                if(elMerged(i) == 0)
                    for q = 1:3
                        j = el2El(i,q);
                        if(j > 0)
                            if(j ~= i && elMerged(i) == 0 && elMerged(j) == 0)
                                for r = 1:3
                                    k = el2El(j,r);
                                    if(k > 0)
                                        if(k ~= i && k ~= j && elMerged(i) == 0 && elMerged(j) == 0 && elMerged(k) == 0)
                                            for s = 1:3
                                                p = el2El(k,s);
                                                if(p > 0)
                                                    if(p ~= i && p ~= j && p ~= k && elMerged(i) == 0 && elMerged(j) == 0 && elMerged(k) == 0 && elMerged(p) == 0)
                                                        allNds = [obj.triElements(i,:),obj.triElements(j,:),obj.triElements(k,:),obj.triElements(p,:)];
                                                        allNds = sort(allNds);
                                                        commonNds = [];
                                                        for m = 1:11
                                                            if(allNds(m) == allNds(m+1))
                                                                commonNds = [commonNds,allNds(m)];
                                                            end
                                                        end
                                                        if(length(commonNds) == 7)
                                                            common2 = [];
                                                            for m = 1:6
                                                                if(commonNds(m) == commonNds(m+1))
                                                                    common2 = [common2,commonNds(m)];
                                                                end
                                                            end
                                                            if(length(common2) == 2)
                                                                if(common2(1) == common2(2))
                                                                    center = common2(1);
                                                                end
                                                            else
                                                                center = 0;
                                                            end
                                                            if(center ~= 0)
                                                                elMerged(i) = 1;
                                                                elMerged(j) = 1;
                                                                elMerged(k) = 1;
                                                                elMerged(p) = 1;
                                                                nodeElim(center) = 1;
                                                                threeEls = sort([obj.triElements(i,:),obj.triElements(j,:),obj.triElements(k,:)]);
                                                                nds1n2 = [];
                                                                for m = 1:8
                                                                    if(threeEls(m) ~= center && threeEls(m) == threeEls(m+1))
                                                                        nds1n2 = [nds1n2,threeEls(m)];
                                                                    end
                                                                end
                                                                nds3n4 = [];
                                                                for m = 1:7
                                                                    mNd = commonNds(m);
                                                                    if(mNd ~= center && mNd ~= nds1n2(1) && mNd ~= nds1n2(2))
                                                                        nds3n4 = [nds3n4,mNd];
                                                                    end
                                                                end
                                                                n1 = nds1n2(1);
                                                                n2 = nds1n2(2);
                                                                n3 = nds3n4(1);
                                                                n4 = nds3n4(2);
                                                                v1 = obj.nodes(n2,:) - obj.nodes(n1,:);
                                                                v2 = obj.nodes(n3,:) - obj.nodes(n2,:);
                                                                v3 = obj.nodes(n4,:) - obj.nodes(n3,:);
                                                                cp1 = v1(1)*v2(2) - v1(2)*v2(1);
                                                                cp2 = v2(1)*v3(2) - v2(2)*v3(1);
                                                                if(cp1*cp2 > 0)
                                                                    newEl = [n1,n2,n3,n4];
                                                                else
                                                                    newEl = [n1,n2,n4,n3];
                                                                end
                                                                obj.quadElements = [obj.quadElements;newEl];
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            newElMerged = elMerged;
            newNodeElim = nodeElim;
        end
        
        function [newElMerged,newNodeElim,newEl2El,obj] = mergeNestedTriEls(obj,elMerged,el2El,nodeElim)
            numEls = size(obj.triElements,1);
            for i = 1:numEls
                if(elMerged(i) == 0)
                    for p = 1:3
                        j = el2El(i,p);
                        if(j > 0)
                            if(j ~= i && elMerged(i) == 0 && elMerged(j) == 0)
                                for q = 1:3
                                    k = el2El(j,q);
                                    if(k > 0)
                                        if(k ~= i && k ~= j && elMerged(i) == 0 && elMerged(j) == 0 && elMerged(k) == 0)
                                            allNds = [obj.triElements(i,:),obj.triElements(j,:),obj.triElements(k,:)];
                                            allNds = sort(allNds);
                                            commonNds = [];
                                            for m = 1:8
                                                if(allNds(m) == allNds(m+1))
                                                    commonNds = [commonNds,allNds(m)];
                                                end
                                            end
                                            if(length(commonNds) == 5)
                                                center = 0;
                                                for m = 1:4
                                                    if(commonNds(m) == commonNds(m+1))
                                                        center = commonNds(m);
                                                    end
                                                end
                                                if(center ~= 0)
                                                    elMerged(i) = 1;
                                                    elMerged(j) = 1;
                                                    elMerged(k) = 1;
                                                    nodeElim(center) = 1;
                                                    newEl = [];
                                                    for m = 1:5
                                                        if(commonNds(m) ~= center)
                                                            newEl = [newEl,commonNds(m)];
                                                        end
                                                    end
                                                    obj.triElements = [obj.triElements;newEl];
                                                    elMerged = [elMerged;0];
                                                    allNbrs = [el2El(i,:),el2El(j,:),el2El(k,:)];
                                                    e2e = [];
                                                    for m = 1:9
                                                        nb = allNbrs(m);
                                                        if(nb ~= i && nb ~= j && nb ~= k)
                                                            e2e = [e2e,nb];
                                                        end
                                                    end
                                                    ln = length(e2e);
                                                    if(ln ~=3)
                                                        e2e = [e2e,zeros(1,3-ln)];
                                                    end
                                                    el2El = [el2El;e2e];
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            newElMerged = elMerged;
            newNodeElim = nodeElim;
            newEl2El = el2El;
        end
        
        function obj = mergeTriEls(obj,elType)
            elMerged = zeros(size(obj.triElements,1),1);
            nodeElim = zeros(length(obj.nodes),1);
            el2El = zeros(size(obj.triElements,1),3);
            for i = 1:size(obj.edges,1)
                if(obj.edges(i,3) > 0 && obj.edges(i,4) > 0)
                    e1 = obj.edges(i,3);
                    e2 = obj.edges(i,4);
                    inserted = 0;
                    j = 1;
                    while(j <= 3 && inserted == 0)
                        if(el2El(e1,j) == 0)
                            el2El(e1,j) = e2;
                            inserted = 1;
                        elseif(el2El(e1,j) == e2)
                            inserted = 1;
                        end
                        j = j + 1;
                    end
                    inserted = 0;
                    j = 1;
                    while(j <= 3 && inserted == 0)
                        if(el2El(e2,j) == 0)
                            el2El(e2,j) = e1;
                            inserted = 1;
                        elseif(el2El(e2,j) == e1)
                            inserted = 1;
                        end
                        j = j + 1;
                    end
                end
            end
            obj.quadElements = [];

            [elMerged,nodeElim,obj] = obj.mergeNestedQuadEls(elMerged,el2El,nodeElim);
            [elMerged,nodeElim,el2El,obj] = obj.mergeNestedTriEls(elMerged,el2El,nodeElim);
            
            if(contains(elType,'quad'))
                for i = 1:size(obj.triElements,1)
                    if(elMerged(i) == 0)
                        [sN1,ang] = obj.sortNodesByAngle(obj.triElements(i,:));
                        for k = 1:3
                            j = el2El(i,k);
                            if(j > 0)
                                if(i ~= j && elMerged(j) == 0 && elMerged(i) == 0)
                                    [sN2,ang] = obj.sortNodesByAngle(obj.triElements(j,:));
                                    if((sN1(2) == sN2(2) && sN1(3) == sN2(3)) || (sN1(2) == sN2(3) && sN1(3) == sN2(2)))
                                        elConn = [sN1(1:2),sN2(1),sN1(3)];
                                        obj.quadElements = [obj.quadElements;elConn];
                                        elMerged(i) = 1;
                                        elMerged(j) = 1;
                                    end
                                end
                            end
                        end
                    end
                end

                [elMerged,obj] = obj.mergeElsBetweenAngles(elMerged,el2El,0.25*pi,0.75*pi);
                [elMerged,obj] = obj.mergeElsBetweenAngles(elMerged,el2El,0.16666*pi,0.83333*pi);
            end

            newNds = [];
            newLabel = zeros(length(obj.nodes),1);
            lab = 0;
            for i = 1:length(obj.nodes)
                if(nodeElim(i) == 0)
                    lab = lab + 1;
                    newLabel(i) = lab;
                    newNds = [newNds;obj.nodes(i,:)];
                end
            end
            obj.nodes = newNds;
            if(contains(elType,'quad'))
                for i = 1:size(obj.triElements,1)
                    if(elMerged(i) == 0)
                        elConn = [obj.triElements(i,:),0];
                        obj.quadElements = [obj.quadElements;elConn];
                    end
                end
                for i = 1:size(obj.quadElements,1)
                    for j = 1:4
                        orig = obj.quadElements(i,j);
                        if(orig ~= 0)
                            obj.quadElements(i,j) = newLabel(orig);
                        end
                    end
                end
            else
                newTriEls = [];
                for i = 1:size(obj.triElements,1)
                    if(elMerged(i) == 0)
                        newTriEls = [newTriEls;obj.triElements(i,:)];
                    end
                end
                for i = 1:size(obj.quadElements,1)
                    nds = obj.quadElements(i,:);
                    triels = [nds(1),nds(2),nds(4);nds(3),nds(4),nds(2)];
                    newTriEls = [newTriEls;triels];
                end
                obj.triElements = newTriEls;
                for i = 1:size(obj.triElements,1)
                    for j = 1:3
                        orig = obj.triElements(i,j);
                        if(orig ~= 0)
                            obj.triElements(i,j) = newLabel(orig);
                        end
                    end
                end
            end
        end
        
        function obj = prepareForMesh(obj)
            numNds = length(obj.nodes);
            numEdges = length(obj.edges);
            obj.triElements = [];
            obj.quadElements = [];
            minSpacing = 10*max(abs(obj.nodes),[],'all');
            maxSpacing = 0;
            avgSpacing = 0;
            for i = 1:numEdges
                n1 = obj.edges(i,1);
                n2 = obj.edges(i,2);
                vec = obj.nodes(n2,:) - obj.nodes(n1,:);
                mag = sqrt(vec*vec');
                avgSpacing = avgSpacing + mag;
                if(mag < minSpacing)
                    minSpacing = mag;
                end
                if(mag > maxSpacing)
                    maxSpacing = mag;
                end
            end
            obj.maxEdgeLen = maxSpacing;
            obj.avgProjLen = 0.866025403784439*(avgSpacing/size(obj.edges,1));
            obj.maxElSize = 0.9*maxSpacing;
            obj.minElSize = 0.8*minSpacing;
            
            minX = min(obj.nodes(:,1),[],'all') - maxSpacing;
            maxX = max(obj.nodes(:,1),[],'all') + maxSpacing;
            minY = min(obj.nodes(:,2),[],'all') - maxSpacing;
            maxY = max(obj.nodes(:,2),[],'all') + maxSpacing;
            avgSpacing = 0.5*(minSpacing+maxSpacing);
            obj.nodeGL = spatialGridList2D(minX,maxX,minY,maxY,avgSpacing,avgSpacing);
            for i = 1:numNds
                obj.nodeGL = obj.nodeGL.addEntry(i,obj.nodes(i,:));
            end
            obj.edgeGL = spatialGridList2D(minX,maxX,minY,maxY,avgSpacing,avgSpacing);
            for i = 1:numEdges
                n1 = obj.edges(i,1);
                n2 = obj.edges(i,2);
                midPt = 0.5*(obj.nodes(n1,:) + obj.nodes(n2,:));
                obj.edgeGL = obj.edgeGL.addEntry(i,midPt);
            end
            obj.triElGL = spatialGridList2D(minX,maxX,minY,maxY,avgSpacing,avgSpacing);
            obj = obj.getBoundaryEdgeNormals(minX,maxX,minY,maxY,sqrt(0.05)*minSpacing);
%             gridSize = minSpacing/4;
%             minX = min(obj.nodes(:,1),[],'all') - 10*gridSize;
%             maxX = max(obj.nodes(:,1),[],'all') + 10*gridSize;
%             minY = min(obj.nodes(:,2),[],'all') - 10*gridSize;
%             maxY = max(obj.nodes(:,2),[],'all') + 10*gridSize;
%             obj.grdStat = gridStatus2D(minX,maxX,minY,maxY,gridSize);
%             obj.grdStat = obj.grdStat.getInteriorGrid(obj.nodes,obj.edges);
%             obj.edges = obj.grdStat.getBoundEdgeNormals(obj.edges,obj.nodes);
        end
        
        function [inEl] = ptInEl(obj,pt,elConn,nodes,marginFactor)
            inEl = 1;
            n1 = elConn(1);
            n2 = elConn(2);
            n3 = elConn(3);
            elNd = [nodes(n1,:);nodes(n2,:);nodes(n3,:)];
            centroid = 0.33333333*(elNd(1,:) + elNd(2,:) + elNd(3,:));
            cenMat = ones(3,2)*[centroid(1),0;0,centroid(2)];
            elNd = elNd - cenMat;
            elNd = marginFactor*elNd;
            elNd = elNd + cenMat;
            v1 = elNd(2,:) - elNd(1,:);
            v2 = elNd(3,:) - elNd(1,:);
            v3 = pt - elNd(1,:);
            cp1 = v1(1)*v2(2) - v1(2)*v2(1);
            cp2 = v1(1)*v3(2) - v1(2)*v3(1);
            if(cp1*cp2 <= 0)
                inEl = 0;
            end
            if(inEl == 1)
                v1 = elNd(3,:) - elNd(2,:);
                v2 = elNd(1,:) - elNd(2,:);
                v3 = pt - elNd(2,:);
                cp1 = v1(1)*v2(2) - v1(2)*v2(1);
                cp2 = v1(1)*v3(2) - v1(2)*v3(1);
                if(cp1*cp2 <= 0)
                    inEl = 0;
                end
            end
            if(inEl == 1)
                v1 = elNd(1,:) - elNd(3,:);
                v2 = elNd(2,:) - elNd(3,:);
                v3 = pt - elNd(3,:);
                cp1 = v1(1)*v2(2) - v1(2)*v2(1);
                cp2 = v1(1)*v3(2) - v1(2)*v3(1);
                if(cp1*cp2 <= 0)
                    inEl = 0;
                end
            end
        end
        
        function obj = skewNodes(obj)
            skewMat = [1,tan(pi/12);tan(pi/12),1];
            obj.nodes = obj.nodes*skewMat;
        end

        function obj = sortIncompleteEdges(obj)
            completeEdges = [];
            completeLabels = [];
            incompleteEdges = [];
            incompleteLabels = [];
            for i = 1:length(obj.edges)
                if(obj.edges(i,4) == 0)
                    incompleteEdges = [incompleteEdges;obj.edges(i,:)];
                    incompleteLabels = [incompleteLabels;i];
                else
                    completeEdges = [completeEdges;obj.edges(i,:)];
                    completeLabels = [completeLabels;i];
                end
            end
            numInc = size(incompleteEdges,1);
            eLen = [];
            for i = 1:numInc
                n1 = incompleteEdges(i,1);
                n2 = incompleteEdges(i,2);
                vec = obj.nodes(n1,1:2) - obj.nodes(n2,1:2);
                mag = sqrt(vec*vec');
                eLen = [eLen;mag];
            end
            for i = 1:numInc
                for j = 1:numInc-1
                    if(eLen(j+1) < eLen(j))
                        edgeSwap = incompleteEdges(j,:);
                        incompleteEdges(j,:) = incompleteEdges(j+1,:);
                        incompleteEdges(j+1,:) = edgeSwap;
                        swap = eLen(j);
                        eLen(j) = eLen(j+1);
                        eLen(j+1) = swap;
                        labSwap = incompleteLabels(j);
                        incompleteLabels(j) = incompleteLabels(j+1);
                        incompleteLabels(j+1) = labSwap;
                    end
                end
            end
            obj.edges = [incompleteEdges;completeEdges];
            sortedLabels = [incompleteLabels;completeLabels];
            obj.edgeGL = obj.edgeGL.reOrderLabels(sortedLabels);
        end

        function [sortedNodes,angles] = sortNodesByAngle(obj,labels)
            l1 = labels(1);
            l2 = labels(2);
            l3 = labels(3);
            angles = [0,0,0];
            v1 = obj.nodes(l2,:) - obj.nodes(l1,:);
            mag1 = sqrt(v1*v1');
            n1 = (1/mag1)*v1;
            v2 = obj.nodes(l3,:) - obj.nodes(l1,:);
            mag2 = sqrt(v2*v2');
            n2 = (1/mag2)*v2;
            dp = n1*n2';
            angles(1) = acos(dp);
            v1 = obj.nodes(l1,:) - obj.nodes(l2,:);
            mag1 = sqrt(v1*v1');
            n1 = (1/mag1)*v1;
            v2 = obj.nodes(l3,:) - obj.nodes(l2,:);
            mag2 = sqrt(v2*v2');
            n2 = (1/mag2)*v2;
            dp = n1*n2';
            angles(2) = acos(dp);
            v1 = obj.nodes(l1,:) - obj.nodes(l3,:);
            mag1 = sqrt(v1*v1');
            n1 = (1/mag1)*v1;
            v2 = obj.nodes(l2,:) - obj.nodes(l3,:);
            mag2 = sqrt(v2*v2');
            n2 = (1/mag2)*v2;
            dp = n1*n2';
            angles(3) = acos(dp);
            sortedNodes = labels;
            for i = 1:3
                for j = 1:2
                    if(angles(j+1) > angles(j))
                        swap = sortedNodes(j);
                        sortedNodes(j) = sortedNodes(j+1);
                        sortedNodes(j+1) = swap;
                        swap = angles(j);
                        angles(j) = angles(j+1);
                        angles(j+1) = swap;
                    end
                end
            end
        end
        
        function obj = uniformBoundarySpacing(obj,maxIt)
            numNds = length(obj.nodes);
            numEdges = length(obj.edges);
            eqnFact = 10;

            cols = 2*numNds;
            rows = 2*cols;
            eqnMat = zeros(rows,cols);
            for i = 1:numNds
                j = 2*i;
                eqnMat(j-1,j-1) = 1;
                eqnMat(j,j) = 1;
            end
            for i = 1:numEdges
                n1 = obj.edges(i,1);
                n2 = obj.edges(i,2);
                j = cols + 2*i;
                k11 = 2*n1 - 1;
                k12 = 2*n1;
                k21 = 2*n2 - 1;
                k22 = 2*n2;
                eqnMat(j-1,k21) = eqnFact;
                eqnMat(j-1,k11) = -eqnFact;
                eqnMat(j,k22) = eqnFact;
                eqnMat(j,k12) = -eqnFact;
            end

            [Q,R] = qr(eqnMat,0);

            for it = 1:maxIt
                edgeDir = obj.getBoundaryData();

                rhs = zeros(rows,1);
                for i = 1:numNds
                    j = 2*i;
                    rhs(j-1) = obj.nodes(i,1);
                    rhs(j) = obj.nodes(i,2);
                end
                for i = 1:numEdges
                    j = cols + 2*i;
                    rhs(j-1) = eqnFact*edgeDir(i,1);
                    rhs(j) = eqnFact*edgeDir(i,2);
                end

                solnVec = R\Q'*rhs;

                obj.nodes = [];
                for i = 1:numNds
                    j = i*2;
                    nd = solnVec(j-1:j)';
                    obj.nodes = [obj.nodes;nd];
                end
            end
        end
        
        function obj = unSkewNodes(obj)
            skewMat = [1,tan(pi/12);tan(pi/12),1];
            invSkew = skewMat\[1,0;0,1];
            obj.nodes = obj.nodes*invSkew;
        end
        
    end
end

