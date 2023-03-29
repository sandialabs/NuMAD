import numpy as np
from pynumad.utils.interpolation import interpolator_wrap


class NuMesh2D():
    """object representing a 2D mesh of triangular or quadrilateral elements
    """

    def __init__(self,boundaryNodes,boundaryEdges): 
        self.nodes = np.array([])
        self.edges = np.array([])
        self.triElements = np.array([])
        self.quadElements = np.array([])
        self.nodeGL = None
        self.edgeGL = None
        self.triElGL = None
        self.maxEdgeLen = 1
        self.avgProjLen = 1
        self.maxElSize = 1
        self.minElSize = 1
        self.nodes = boundaryNodes
        numNds = boundaryNodes.shape[0]
        if not boundaryEdges.any():
            self.edges = np.zeros((numNds,6))
            self.edges[:,0] = np.array([range(numNds)])
            self.edges[:,1] = np.array([range(1,numNds),1])
            self.edges[:,2] = - 1
        else:
            numEdges = boundaryEdges.shape[0]
            self.edges = np.hstack([boundaryEdges,- np.ones((numEdges,1)),np.zeros((numEdges,3))])
        
        return
        
        
    def adoptAnyNode(self, currentEdge, ndPt, srchRad):
        """Object data modified: self.triElements, self.triElGL, self.edges

        Parameters
        ----------
        currentEdge
        ndPt
        srchRad

        Returns
        -------
        nodeFound
        edgesAdded
        """
        nodeFound = 0
        edgesAdded = 0
        n1 = self.edges[currentEdge,0]
        n2 = self.edges[currentEdge,1]
        edgeMpt = 0.5 * (self.nodes[n1,:] + self.nodes[n2,:])
        nearNds = self.nodeGL.findInRadius(ndPt,srchRad)
        k = 1
        while (nodeFound == 0 and k <= len(nearNds)):
            j = nearNds[k]
            if (j != n1 and j != n2):
                vec = self.nodes[j,:] - ndPt
                dist = np.linalg.norm(vec)
                mpVec = self.nodes[j,:] - edgeMpt
                dp = self.edges[currentEdge,4:6] * mpVec.T
                if (dist < srchRad and dp > 0):
                    elNum = self.findElement(np.array([n1,n2,j]))
                    if (elNum == 0):
                        violation = self.checkViolations(np.array([n1,n2,j]),[],1.3)
                        if (violation == 0):
                            numAdded,self = self.completeEdges(n1,n2,j,self.triElements.shape[0] + 1)
                            edgesAdded = edgesAdded + numAdded
                            newEl = np.array([n1,n2,j])
                            self.triElements = np.array([[self.triElements],[newEl]])
                            midPt = 0.333333 * (self.nodes[n1,:] + self.nodes[n2,:] + self.nodes[j,:])
                            self.triElGL = self.triElGL.addEntry(self.triElements.shape[0]-1,midPt)
                            self.edges[currentEdge,3] = self.triElements.shape[0]-1
                            nodeFound = 1
                    else:
                        self.edges[currentEdge,3] = elNum
                        nodeFound = 1
            k = k+1
        return nodeFound,edgesAdded
        
        
    def adoptConnectedNode(self, currentEdge, ndPt, srchRad):
        """Object data modified: self.triElements, self.triElGL, self.edges,

        Parameters
        ----------
        elType
        equalizeSpacing

        Returns
        -------
        nodeFound
        edgesAdded
        """
        nodeFound = 0
        edgesAdded = 0
        n1 = self.edges[currentEdge,0].astype(int)
        n2 = self.edges[currentEdge,1].astype(int)
        edgeMpt = 0.5 * (self.nodes[n1,:] + self.nodes[n2,:])
        nearEdges = self.edgeGL.findInRadius(edgeMpt,1.1 * self.maxEdgeLen)
        k = 1
        while (nodeFound == 0 and k <= len(nearEdges)):
            j = nearEdges[k]
            if (currentEdge != j):
                n1j = self.edges[j,0].astype(int)
                n2j = self.edges[j,0].astype(int)
                if (n1j == n1 or n1j == n2 or n2j == n1 or n2j == n2):
                    if (n1j == n1 or n1j == n2):
                        n3 = n2j
                    else:
                        n3 = n1j
                    vec = self.nodes[n3,0:2] - ndPt
                    dist = np.linalg.norm(vec)
                    mpVec = self.nodes[n3,:] - edgeMpt
                    dp = self.edges[currentEdge,4:6] @ mpVec.T
                    if (dist < srchRad and dp > 0):
                        elNum = self.findElement(np.array([n1,n2,n3]))
                        if (elNum == 0):
                            violation = self.checkViolations(np.array([n1,n2,n3]),[],1.3)
                            if (violation == 0):
                                numAdded,self = self.completeEdges(n1,n2,n3,self.triElements.shape[0] + 1)
                                edgesAdded = edgesAdded + numAdded
                                newEl = np.array([n1,n2,n3])
                                self.triElements = np.array([[self.triElements],[newEl]])
                                midPt = 0.333333 * (self.nodes[n1,:] + self.nodes[n2,:] + self.nodes[n3,:])
                                self.triElGL = self.triElGL.addEntry(self.triElements.shape[0]-1,midPt)
                                self.edges[currentEdge,3] = self.triElements.shape[0]-1
                                nodeFound = 1
                        else:
                            self.edges[currentEdge,3] = elNum
                            nodeFound = 1
            k = k+1
        return nodeFound,edgesAdded
        
        
    def checkViolations(self, ndLab, ndCrd, marginFact):
        """Object data modified: none
        Parameters
        ----------
        ndLab
        ndCrd
        marginFact

        Returns
        -------
        violation : bool
        """ 
        violation = 0
        elCrd = []
        for i in range(len(ndLab)):
            if (ndLab[i] > 0):
                elCrd = np.array([[elCrd],[self.nodes[ndLab[i],:]]])
            else:
                j = - ndLab[i]
                elCrd = np.array([[elCrd],[ndCrd[j,:]]])
        
        midPt = 0.333333 * (elCrd[0,:] + elCrd[1,:] + elCrd[2,:])
        nearEls = self.triElGL.findInRadius(midPt,2.1 * self.maxEdgeLen)
        for m in range(len(nearEls)):
            nm = self.triElements[nearEls[m],:]
            thisElnds = np.array([[self.nodes[nm[0],:]],[self.nodes[nm[1],:]],[self.nodes[nm[3],:]]])
            overlap = self.elementsOverlap(np.array([[1,2,3],[4,5,6]]),np.array([[elCrd],[thisElnds]]),1)
            if (overlap == 1):
                violation = 1
            for i in range(3):
                if (ndLab[i] < 0):
                    inEl = self.ptInEl(elCrd[i,:],np.array([1,2,3]),thisElnds,marginFact)
                    if (inEl == 1):
                        violation = 1
        
        nearNds = self.nodeGL.findInRadius(midPt,1.1 * self.maxEdgeLen)
        for m in range(len(nearNds)):
            nd = nearNds[m]
            if (nd != ndLab[0] and nd != ndLab[1] and nd != ndLab[2]):
                pt = self.nodes[nd,:]
                inEl = self.ptInEl(pt,np.array([1,2,3]),elCrd,marginFact)
                if (inEl == 1):
                    violation = 1
        return violation
        
        
    def completeEdges(self, n1, n2, n3, elNum):
        """Object data modified: self.edges, self.edgeGL, self.maxEdgeLen
        Parameters
        ----------
        n1
        n2
        n3
        elNum

        Returns
        -------
        numAdded
        """ 
        numAdded = 0
        midPt = 0.333333 * (self.nodes[n1,:] + self.nodes[n2,:] + self.nodes[n3,:])
        nearEdges = self.edgeGL.findInRadius(midPt,1.1 * self.maxEdgeLen)
        e2Found = 0
        j = 1
        while (e2Found == 0 and j <= len(nearEdges)):
    
            k = nearEdges[j]
            if ((self.edges[k,1] == n1 and self.edges[k,2] == n3) or (self.edges[k,1] == n3 and self.edges[k,2] == n1)):
                e2Found = k
            j = j+1
    
        if (e2Found == 0):
            newEdge = np.array([n1,n3,elNum,0,0,0])
            self.edges = np.vstack([self.edges,newEdge])
            numAdded = numAdded + 1
            midPt = 0.5 * (self.nodes[n1,:] + self.nodes[n3,:])
            self.edgeGL = self.edgeGL.addEntry(self.edges.shape[0]-1,midPt)
            vec = self.nodes[n1,:] - self.nodes[n3,:]
            mag = np.linalg.norm(vec)
            if (mag > self.maxEdgeLen):
                self.maxEdgeLen = mag
        else:
            self.edges[e2Found,3] = elNum
        
        e3Found = 0
        j = 1
        while (e3Found == 0 and j <= len(nearEdges)):
    
            k = nearEdges[j]
            if ((self.edges[k,0] == n2 and self.edges[k,1] == n3) or 
                    (self.edges[k,0] == n3 and self.edges[k,1] == n2)):
                e3Found = k
            j = j+1
        
        if (e3Found == 0):
            newEdge = np.array([n2,n3,elNum,0,0,0])
            self.edges = np.vstack([self.edges,newEdge])
            numAdded = numAdded + 1
            midPt = 0.5 * (self.nodes[n2,:] + self.nodes[n3,:])
            self.edgeGL = self.edgeGL.addEntry(self.edges.shape[0]-1,midPt)
            vec = self.nodes[n2,:] - self.nodes[n3,:]
            mag = np.linalg.norm(vec)
            if (mag > self.maxEdgeLen):
                self.maxEdgeLen = mag
        else:
            self.edges[e3Found,3] = elNum
        
        return numAdded
        
        
    def createNode(self, currentEdge, ndPt):
        """Object data modified: self.nodes, self.nodeGL, self.triElements,
        self.triGL
        Parameters
        ----------

        Returns
        -------
        nodeCreated
        edgesAdded
        """ 
        nodeCreated = 0
        edgesAdded = 0
        n1 = self.edges[currentEdge,0]
        n2 = self.edges[currentEdge,1]
        violation = self.checkViolations(np.array([n1,n2,-1]),ndPt,1.3)
        if (violation == 0):
            newNd = ndPt
            self.nodes = np.array([[self.nodes],[newNd]])
            self.nodeGL = self.nodeGL.addEntry(self.nodes.shape[0]-1,newNd)
            n3 = self.nodes.shape[0]
            newEl = np.array([n1,n2,n3])
            self.triElements = np.array([[self.triElements],[newEl]])
            midPt = 0.333333 * (self.nodes[n1,:] + self.nodes[n2,:] + self.nodes[n3,:])
            self.triElGL = self.triElGL.addEntry(self.triElements.shape[0]-1,midPt)
            numEls = self.triElements.shape[0]
            numAdded,self = self.completeEdges(n1,n2,n3,numEls)
            edgesAdded = edgesAdded + numAdded
            nodeCreated = 1
        
        return nodeCreated,edgesAdded
        
        
    def createPlanarMesh(self, elType, equalizeSpacing):
        """Object data modified: none
        Parameters
        ----------
        elType
        equalizeSpacing

        Returns
        -------
        nodes
        elements
        """ 
        if ('quad' in elType):
            self = self.skewNodes()
        
        boundaryNodes = self.nodes
        if (equalizeSpacing):
            maxIt = np.ceil(0.5 * self.nodes.shape[0]).astype(int)
            self = self.uniformBoundarySpacing(maxIt)
        
        #             plot2DMesh(self.nodes,[]);
    #             keyboard
        self.createUnstructTriMesh()
        #             plot2DMesh(self.nodes,self.triElements);
    #             keyboard
        self.distributeNodes(boundaryNodes)
        #             plot2DMesh(self.nodes,self.triElements);
    #             keyboard
        if ('quad' in elType):
            self.unSkewNodes()
            #                 plot2DMesh(self.nodes,self.triElements);
    #                 keyboard
        
        self.mergeTriEls(elType)
        nodes = self.nodes
        if ('quad' in elType):
            elements = self.quadElements
        else:
            elements = self.triElements
        
        #             plot2DMesh(nodes,elements);
    #             keyboard;
        return nodes,elements
        
        
    def createSweptMesh(self, sweepMethod, direction, sweepDistance, sweepElements, followNormal, destNodes): 
        """Object data modified: self.quadElements, self.nodes, self.quadElements
        Parameters
        ----------

        Returns
        -------
        nodes
        elements
        """ 
        numNds = self.nodes.shape[0]
        numEdges = self.edges.shape[0]
        self.quadElements = np.array([])
        if (followNormal == 1):
            methString = 'in_direction, to_point, from_point'
            if (sweepMethod in methString):
                pass
        else:
            if 'in_direction' == sweepMethod:
                dirMag = np.linalg.norm(direction)
                unitDir = (1 / dirMag) * direction
                ndDir = (sweepDistance / sweepElements) * np.matmul(np.ones((numNds,2)), np.array([[unitDir[0],0],[0,unitDir[1]]]))
        
        if ('revolve' in sweepMethod):
            pass
        else:
            for i in range(1,sweepElements+1):
                ndRow = self.nodes[:numNds,:] + i * ndDir
                self.nodes = np.concatenate([self.nodes,ndRow])
                for j in range(numEdges):
                    n1 = (i - 1) * numNds + self.edges[j,0]
                    n2 = (i - 1) * numNds + self.edges[j,1]
                    n3 = i * numNds + self.edges[j,1]
                    n4 = i * numNds + self.edges[j,0]
                    try:
                        self.quadElements = np.vstack((self.quadElements,[n1,n2,n3,n4]))
                    except ValueError:
                        self.quadElements = [n1,n2,n3,n4]
        
        nodes = self.nodes
        elements = self.quadElements.astype(int)
        return nodes,elements
        
        
    def createUnstructTriMesh(self):
        """Object data modified: self.triElements, self.triElGL,
        self.edges, self.edgeGL, self.maxEdgeLen

        Parameters
        ----------

        Returns
        -------
        self
        """ 
        self.prepareForMesh()
        numEdges = self.edges.shape[0]
        for i in range(numEdges):
            n1i = self.edges[i,0].astype(int)
            n2i = self.edges[i,1].astype(int)
            midPt = 0.5 * (self.nodes[n1i,:] + self.nodes[n2i,:])
            nearEdges = self.edgeGL.findInRadius(midPt,1.1 * self.maxEdgeLen)
            for m in range(len(nearEdges)):
                j = nearEdges[m]
                if (j != i):
                    n1j = self.edges[j,0].astype(int)
                    n2j = self.edges[j,1].astype(int)
                    allNds = np.array([n1i,n2i,n1j,n2j])
                    allNds = np.sort(allNds)
                    n1 = 0
                    for k in range(3):
                        if (allNds[k+1] == allNds[k]):
                            n1 = allNds[k]
                            allNds[k] = 0
                            allNds[k+1] = 0
                    allNds = np.sort(allNds)
                    n2 = allNds[2]
                    n3 = allNds[3]
                    if (n1 != 0):
                        v1 = self.nodes[n2,:] - self.nodes[n1,:]
                        mag1 = np.linalg.norm(v1)
                        v2 = self.nodes[n3,:] - self.nodes[n1,:]
                        mag2 = np.linalg.norm(v2)
                        dp = (v1 @ v2.T) / (mag1 * mag2)
                        if (dp > 0.4):
                            v3 = v1 + v2
                            dp = self.edges[i,4:6] @ v3
                            if (dp > 0):
                                newEl = np.array([n1,n2,n3])
                                eNum = self.findElement(newEl)
                                if (eNum == 0):
                                    self.triElements = np.vstack([self.triElements,newEl]) if self.triElements.size else newEl.reshape(1,-1)
                                    midPt = 0.3333333 * (self.nodes[n1,:] + self.nodes[n2,:] + self.nodes[n3,:])
                                    self.triElGL = self.triElGL.addEntry(self.triElements.shape[0]-1,midPt)
                                    eNum = self.triElements.shape[0] - 1
                                    self.edges[i,3] = eNum
                                    self.edges[j,3] = eNum
                                    newEdge = np.array([n2,n3,eNum,0,0,0])
                                    self.edges = np.vstack([self.edges,newEdge])
                                    midPt = 0.5 * (self.nodes[n2,:] + self.nodes[n3,:])
                                    self.edgeGL = self.edgeGL.addEntry(self.edges.shape[0]-1,midPt)
                                    vec = self.nodes[n2,:] - self.nodes[n3,:]
                                    mag = np.linalg.norm(vec)
                                    if (mag > self.maxEdgeLen):
                                        self.maxEdgeLen = mag
        
        edgesAdded = self.edges.shape[0]
        maxIt = edgesAdded
        it = 0
        while (edgesAdded > 0 and it < maxIt):
    
            edgesAdded = 0
            self.sortIncompleteEdges()
            numEdges = self.edges.shape[0]
            #                 plot2DMesh(self.nodes,self.triElements);
    #                 keyboard
            for i in range(numEdges):
                if (self.edges[i,2] > 0):
                    n1 = self.edges[i,0].astype(int)
                    n2 = self.edges[i,1].astype(int)
                    midpt = 0.5 * (self.nodes[n2,:] + self.nodes[n1,:])
                    eVec = self.nodes[n2,:] - self.nodes[n1,:]
                    eNorm = np.array([- eVec[1],eVec[0]])
                    el = self.edges[i,2].astype(int)
                    for j in range(3):
                        if (self.triElements[el,j] != n1 and self.triElements[el,j] != n2):
                            n3 = self.triElements[el,j]
                            vin = midpt - self.nodes[n3,:]
                            dp = eNorm @ vin
                            if (dp < 0):
                                eNorm = - eNorm
                            mag = np.linalg.norm(eNorm)
                            self.edges[i,4:6] = (1 / mag) * eNorm
            for i in range(numEdges):
                if (self.edges[i,3] == 0):
                    n1 = self.edges[i,0].astype(int)
                    n2 = self.edges[i,1].astype(int)
                    midpt = 0.5 * (self.nodes[n1,0:2] + self.nodes[n2,1:3])
                    vec = self.nodes[n1,:3] - self.nodes[n2,:3]
                    eLen = np.linalg.norm(vec)
                    unitProj = self.edges[i,4:6]
                    projLen = 0.8 * self.avgProjLen + 0.2 * 0.866025403784439 * eLen
                    srchRad = 0.5 * np.sqrt(projLen ** 2 + eLen ** 2)
                    ndPt = midpt + 0.5 * projLen * unitProj
                    nodeFound,edAdd = self.adoptConnectedNode(i,ndPt,srchRad)
                    edgesAdded = edgesAdded + edAdd
                    if (nodeFound == 0):
                        ndPt = midpt + projLen * unitProj
                        nodeFound,edAdd,self = self.adoptConnectedNode(i,ndPt,srchRad)
                        edgesAdded = edgesAdded + edAdd
                    if (nodeFound == 0):
                        ndPt = midpt + 0.5 * projLen * unitProj
                        nodeFound,edAdd,self = self.adoptAnyNode(i,ndPt,srchRad)
                        edgesAdded = edgesAdded + edAdd
                    if (nodeFound == 0):
                        ndPt = midpt + projLen * unitProj
                        nodeFound,edAdd,self = self.adoptAnyNode(i,ndPt,srchRad)
                        edgesAdded = edgesAdded + edAdd
                    if (nodeFound == 0):
                        ndPt = midpt + projLen * unitProj
                        nodeCreated,edAdd,self = self.createNode(i,ndPt)
                        edgesAdded = edgesAdded + edAdd
                        if (nodeCreated == 0):
                            ndPt = midpt + 0.5 * projLen * unitProj
                            nodeCreated,edAdd,self = self.createNode(i,ndPt)
                            edgesAdded = edgesAdded + edAdd
            #                 tic
    #                 for i = 1:length(self.edges)
    #                     if(self.edges(i,4) == 0)
    #                         n1 = self.edges[i,0];
    #                         n2 = self.edges[i,1];
    #                         midPt = 0.5*(self.nodes[n1,:] + self.nodes[n2,:]);
    #                         nearEdges = self.edgeGL.findInRadius(midPt,1.1*self.maxEdgeLen);
    #                         for p = 1:length(nearEdges)
    #                             j = nearEdges[p];
    #                             if(j ~= i && self.edges(j,4) == 0 && self.edges(i,4) == 0)
    #                                 for q = 1:length(nearEdges)
    #                                     k = nearEdges[q];
    #                                     if(k ~= i && k ~= j && self.edges(i,4) == 0 && self.edges(j,4) == 0 && self.edges(k,4) == 0)
    #                                         allNds = [self.edges(i,1:2),self.edges(j,1:2),self.edges(k,1:2)];
    #                                         allNds = sort(allNds);
    #                                         numDup = 0;
    #                                         for m = 1:5
    #                                             if(allNds(m+1) == allNds[m])
    #                                                 numDup = numDup + 1;
    #                                             end
    #                                         end
    #                                         if(numDup == 3)
    #                                             newEl = [allNds[0],allNds[2],allNds(5)];
    #                                             elNum = self.findElement(newEl);
    #                                             if(elNum == 0)
    #                                                 violation = self.checkViolations(newEl,[],1.3);
    #                                                 if(violation == 0)
    #                                                     self.triElements = [self.triElements;newEl];
    #                                                     elNum = size(self.triElements,1);
    #                                                     midPt = 0.333333*(self.nodes(newEl[0],:) + self.nodes(newEl[1],:) + self.nodes(newEl[2],:));
    #                                                     self.triElGL = self.triElGL.addEntry(elNum,midPt);
    #                                                     self.edges(i,4) = elNum;
    #                                                     self.edges(j,4) = elNum;
    #                                                     self.edges(k,4) = elNum;
    #                                                 end
    #                                             end
    #                                         end
    #                                     end
    #                                 end
    #                             end
    #                         end
    #                     end
    #                 end
    #                 toc
            it = it + 1
    
        
        return self
        
        
    def distributeNodes(self, boundaryNodes):
        """Object data modified: self.nodes

        Parameters
        ----------
        boundaryNodes

        Returns
        -------
        self
        """ 
        numNds = self.nodes.shape[0]
        dim = 2 * numNds
        Dmat = np.zeros((dim))
        numBound = boundaryNodes.shape[0]
        Dmat[np.arange(2 * numBound)] = 100000
        Pmat = 10 * np.ones((dim)) + Dmat
        np.pinv = 1.0 / Pmat
        rhs = np.zeros((dim))
        for i in range(numBound):
            j = 2*(i+1) - 1
            rhs[np.arange(j-1,j+1)] = 100000 * boundaryNodes[i,:].T
        
        numEls = self.triElements.shape[0]
        elWt = np.zeros((numEls,1))
        for i in range(numEls):
            ni = self.triElements[i,:]
            v1 = self.nodes[ni[1],:] - self.nodes[ni[0],:]
            v2 = self.nodes[ni[2],:] - self.nodes[ni[0],:]
            cp = v1[0] * v2[1] - v1[1] * v2[0]
            elWt[i] = np.abs(cp)
        
        avgWt = np.mean(elWt,'all')
        elWt = (1 / avgWt) * elWt
        eMat = np.array([[2,0,- 1,0,- 1,0],[0,2,0,- 1,0,- 1],[- 1,0,2,0,- 1,0],[0,- 1,0,2,0,- 1],[- 1,0,- 1,0,2,0],[0,- 1,0,- 1,0,2]])
        xVec = np.zeros((dim,1))
        gVec = - rhs
        wVec = np.multiply(np.pinv,gVec)
        hVec = - wVec
        res = gVec.T * wVec
        i = 0
        while (res > 1e-12 and i < dim):
    
            zVec = np.zeros((dim,1))
            for j in np.arange(1,self.triElements.shape[0]+1):
                for k in range(3):
                    nd = self.triElements[j,k]
                    for m in range(2):
                        ei = 2 * (k-1) + m
                        gi = 2 * (nd-1) + m
                        for n in range(3):
                            ndj = self.triElements[j,n]
                            for p in range(2):
                                ej = 2 * (n-1) + p
                                gj = 2 * (ndj-1) + p
                                zVec[gi] = zVec(gi) + elWt[j] * eMat(ei,ej) * hVec(gj)
            zVec = zVec + np.multiply(Dmat,hVec)
            alpha = res / (hVec.T * zVec)
            xVec = xVec + alpha * hVec
            gVec = gVec + alpha * zVec
            wVec = np.multiply(np.pinv,gVec)
            rNext = gVec.T * wVec
            beta = rNext / res
            res = rNext
            hVec = - wVec + beta * hVec
            i = i + 1
    
        
        for i in range(numNds):
            j = i * 2
            self.nodes[i,:] = np.transpose(xVec[j-1:j])
        
        return self
        
        
    def elementsOverlap(self, els, ndCrd, marginFact):
        """Object data modified: none
        Parameters
        ----------
        els
        ndCrd
        marginFact

        Returns
        -------
        overlap
        """ 
        Xel1 = np.array([[ndCrd[els[0,0],:]],[ndCrd[els[0,1],:]],[ndCrd[els[0,2],:]]])
        Xel2 = np.array([[ndCrd[els[1,0],:]],[ndCrd[els[1,1],:]],[ndCrd[els[1,2],:]]])
        cent = 0.3333333 * (Xel1[0,:] + Xel1[1,:] + Xel1[2,:])
        centMat = np.ones((3,2)) * np.array([[cent[0],0],[0,cent[1]]])
        Xel1 = marginFact * (Xel1 - centMat) + centMat
        elEdges = np.array([[1,2],[2,3],[3,1]])
        overlap = 0
        for i in range(3):
            n11 = elEdges[i,0]
            n21 = elEdges[i,1]
            x01 = np.transpose(Xel1[n11,:])
            v1 = np.transpose((Xel1[n21,:] - Xel1[n11,:]))
            for j in range(3):
                n12 = elEdges[j,0]
                n22 = elEdges[j,1]
                x02 = np.transpose(Xel2[n12,:])
                v2 = np.transpose((Xel2[n22,:] - Xel2[n12,:]))
                vMat = np.array([v1,- v2])
                mag = np.amax(np.abs(vMat))
                Q,R = np.linalg.qr(vMat,0)
                if (np.sqrt(np.abs(R(1,1) * R(2,2))) > 1e-06 * mag):
                    soln = np.linalg.solve(R,Q.T) * (x02 - x01)
                    if (soln[0] > 1e-06 and soln[0] < 0.999999 and soln[1] > 1e-06 and soln[1] < 0.999999):
                        overlap = 1
        
        return overlap
        
        
    def findElement(self, nds):
        """Object data modified: none
        Parameters
        ----------
        nds

        Returns
        -------
        elNum
        """ 
        nearEls = self.triElGL.findInRadius(
            self.nodes[nds[0],:],
            1.1 * self.maxEdgeLen)
        sortedNds = np.sort(nds)
        elNum = 0
        for j in range(len(nearEls)):
            i = nearEls[j]
            sortedEl = np.sort(self.triElements[i,:])
            if (sortedNds == sortedEl).all():
                elNum = i
        
        return elNum
        
        
    def getBoundaryData(self):
        """Object data modified: none
        Parameters
        ----------

        Returns
        -------
        edgeDir
        """  
        numEdges = self.edges.shape[0]
        avgSpacing = 0
        edgeDir = np.array([])
        for i in range(numEdges):
            n1 = self.edges[i,0].astype(int)
            n2 = self.edges[i,1].astype(int)
            vec = self.nodes[n2,:] - self.nodes[n1,:]
            mag = np.linalg.norm(vec)
            avgSpacing = avgSpacing + mag
            unitE = (1 / mag) * vec
            edgeDir = np.vstack([edgeDir,unitE]) if edgeDir.size else unitE
        
        avgSpacing = avgSpacing / numEdges
        edgeDir = avgSpacing * edgeDir
        return edgeDir
        
        
    def getBoundaryEdgeNormals(self, Xmin, Xmax, Ymin, Ymax, spacing):
        """Object data modified: self.edges
        Parameters
        ----------
        Xmin
        Xmax
        Ymin
        Ymax
        spacing

        Returns
        -------
        self
        """ 
        for i in range(self.edges.shape[0]):
            n1 = self.edges[i,0].astype(int)
            n2 = self.edges[i,1].astype(int)
            v2 = self.nodes[n2,:] - self.nodes[n1,:]
            eNorm = np.array([- v2[1],v2[0]])
            mag = np.linalg.norm(eNorm)
            self.edges[i,4:6] = (1 / mag) * eNorm
        
        for x in np.arange(Xmin,Xmax,spacing):
            v1 = np.array([0,1])
            x01 = np.array([x,Ymin])
            intersects = np.array([], dtype=int)
            nearEdges = self.edgeGL.findInXYMargin(x01,0.6 * self.maxEdgeLen,- 1)
            for j in range(len(nearEdges)):
                i = nearEdges[j]
                n1 = self.edges[i,0].astype(int)
                n2 = self.edges[i,1].astype(int)
                v2 = self.nodes[n2,:].T - self.nodes[n1,:].T
                x02 = self.nodes[n1,:].T
                vMat = np.stack([v1,-v2],axis=1)
                mag = np.amax(np.abs(vMat))
                Q,R = np.linalg.qr(vMat)
                if (np.sqrt(np.abs(R[0,0] * R[1,1])) > 1e-06 * mag):
                    soln = np.linalg.solve(R,Q.T) @ (x02 - x01)
                    if (soln[1] > 0 and soln[1] < 1):
                        solnvec = np.array([i,soln[0]], dtype=int)
                        intersects = np.vstack([intersects,solnvec]) if intersects.size else solnvec
            numInt = intersects.shape[0]
            if (np.abs(np.round(0.5 * numInt) - 0.5 * numInt) < 0.001):
                for i in range(numInt):
                    for j in range(numInt-1):
                        if (intersects[j+1,1] < intersects[j,1]):
                            swap = intersects[j+1,:]
                            intersects[j+1,:] = intersects[j,:]
                            intersects[j,:] = swap
                for i in np.arange(0,numInt,2):
                    ed = intersects[i,0]
                    if (self.edges[ed,5] < 0):
                        self.edges[ed,4:6] = - self.edges[ed,4:6]
                for i in np.arange(0,numInt,2):
                    ed = intersects[i,0]
                    if (self.edges[ed,5] > 0):
                        self.edges[ed,4:6] = - self.edges[ed,4:6]
        
        for y in np.arange(Ymin,Ymax+spacing,spacing):
            v1 = np.array([1,0])
            x01 = np.array([Xmin,y])
            intersects = np.array([],dtype=int)
            nearEdges = self.edgeGL.findInXYMargin(x01,- 1,0.6 * self.maxEdgeLen)
            for j in range(len(nearEdges)):
                i = nearEdges[j]
                n1 = self.edges[i,0].astype(int)
                n2 = self.edges[i,1].astype(int)
                v2 = self.nodes[n2,:].T - self.nodes[n1,:].T
                x02 = self.nodes[n1,:].T
                vMat = np.stack([v1,-v2], axis=1)
                mag = np.amax(np.abs(vMat))
                Q,R = np.linalg.qr(vMat)
                if (np.sqrt(np.abs(R[0,0] * R[1,1])) > 1e-06 * mag):
                    soln = np.linalg.solve(R,Q.T) @ (x02 - x01)
                    if (soln[1] > 0 and soln[1] < 1):
                        solnvec = np.array([i,soln[0]], dtype=int)
                        intersects = np.vstack([intersects,solnvec]) if intersects.size else solnvec
            numInt = intersects.shape[0]
            if (np.abs(np.round(0.5 * numInt) - 0.5 * numInt) < 0.001):
                for i in range(numInt):
                    for j in range(numInt-1):
                        if (intersects[j+1,1] < intersects[j,1]):
                            swap = intersects[j+1,:]
                            intersects[j+1,:] = intersects[j,:]
                            intersects[j,:] = swap
                for i in np.arange(0,numInt,2):
                    ed = intersects[i,0]
                    if (self.edges[ed,5] < 0):
                        self.edges[ed,4:6] = - self.edges[ed,4:6]
                for i in np.arange(0,numInt,2):
                    ed = intersects[i,0]
                    if (self.edges[ed,5] > 0):
                        self.edges[ed,4:6] = - self.edges[ed,4:6]
        
        return self
        
        
    def mergeElsBetweenAngles(self, elMerged, el2El, minAngle, maxAngle):
        """Object data modified: self.quadElements
        Parameters
        ----------
        elMerged
        el2El
        minAngle
        maxAngle

        Returns
        -------
        newElMerged
        """ 
        for i in range(self.triElements.shape[0]):
            if (elMerged[i] == 0):
                for p in range(3):
                    j = el2El[i,p]
                    if (j > 0):
                        if (elMerged[i] == 0 and elMerged[j] == 0 and i != j):
                            allNds = np.array([self.triElements[i,:],self.triElements[j,:]])
                            allNds = np.sort(allNds)
                            commonNds = []
                            for k in range(5):
                                if (allNds[k] == allNds[k+1]):
                                    commonNds = np.array([commonNds,allNds[k]])
                                    allNds[k] = 0
                                    allNds[k+1] = 0
                            if (len(commonNds) == 2):
                                allNds = np.sort(allNds)
                                unCommon = allNds[4:6]
                                quadNds = np.array([unCommon[0],commonNds[0],unCommon[1],commonNds[1]])
                                next = np.array([quadNds[1:5],quadNds[0]])
                                last = np.array([quadNds[3],quadNds[:4]])
                                mergeOk = 1
                                for k in range(4):
                                    v1 = self.nodes[next[k],:] - self.nodes[quadNds[k],:]
                                    mag1 = np.sqrt(v1 * v1.T)
                                    v2 = self.nodes[last[k],:] - self.nodes[quadNds[k],:]
                                    mag2 = np.sqrt(v2 * v2.T)
                                    dp = v1 * v2.T
                                    theta = np.arccos(dp / (mag1 * mag2))
                                    if (theta < minAngle or theta > maxAngle):
                                        mergeOk = 0
                                if (mergeOk == 1):
                                    self.quadElements = np.array([[self.quadElements],[quadNds]])
                                    elMerged[i] = 1
                                    elMerged[j] = 1
        
        newElMerged = elMerged
        return newElMerged
        
        
    def mergeNestedQuadEls(self, elMerged, el2El, nodeElim):
        """Object data modified: self.quadElements
        
        Parameters
        ----------
        elMerged
        el2El
        nodeElim

        Returns
        -------
        newElMerged
        newNodeElim
        """ 
        numEls = self.triElements.shape[0]
        for i in range(numEls):
            if (elMerged[i] == 0):
                for q in range(3):
                    j = el2El(i,q)
                    if (j > 0):
                        if (j != i and elMerged[i] == 0 and elMerged[j] == 0):
                            for r in range(3):
                                k = el2El(j,r)
                                if (k > 0):
                                    if (k != i and k != j and elMerged[i] == 0 and elMerged[j] == 0 and elMerged[k] == 0):
                                        for s in range(3):
                                            p = el2El(k,s)
                                            if (p > 0):
                                                if (p != i and p != j and p != k and elMerged[i] == 0 and elMerged[j] == 0 and elMerged[k] == 0 and elMerged[p] == 0):
                                                    allNds = np.array([self.triElements[i,:],self.triElements[j,:],self.triElements[k,:],self.triElements[p,:]])
                                                    allNds = np.sort(allNds)
                                                    commonNds = []
                                                    for m in range(11):
                                                        if (allNds[m] == allNds(m + 1)):
                                                            commonNds = np.array([commonNds,allNds[m]])
                                                    if (len(commonNds) == 7):
                                                        common2 = []
                                                        for m in range(6):
                                                            if (commonNds[m] == commonNds(m + 1)):
                                                                common2 = np.array([common2,commonNds[m]])
                                                        if (len(common2) == 2):
                                                            if (common2[0] == common2[1]):
                                                                center = common2[0]
                                                        else:
                                                            center = 0
                                                        if (center != 0):
                                                            elMerged[i] = 1
                                                            elMerged[j] = 1
                                                            elMerged[k] = 1
                                                            elMerged[p] = 1
                                                            nodeElim[center] = 1
                                                            threeEls = np.sort(np.array([self.triElements[i,:],self.triElements[j,:],self.triElements[k,:]]))
                                                            nds1n2 = []
                                                            for m in range(8):
                                                                if (threeEls[m] != center and threeEls[m] == threeEls(m + 1)):
                                                                    nds1n2 = np.array([nds1n2,threeEls[m]])
                                                            nds3n4 = []
                                                            for m in range(7):
                                                                mNd = commonNds[m]
                                                                if (mNd != center and mNd != nds1n2[0] and mNd != nds1n2[1]):
                                                                    nds3n4 = np.array([nds3n4,mNd])
                                                            n1 = nds1n2[0]
                                                            n2 = nds1n2[1]
                                                            n3 = nds3n4[0]
                                                            n4 = nds3n4[1]
                                                            v1 = self.nodes[n2,:] - self.nodes[n1,:]
                                                            v2 = self.nodes[n3,:] - self.nodes[n2,:]
                                                            v3 = self.nodes[n4,:] - self.nodes[n3,:]
                                                            cp1 = v1[0] * v2[1] - v1[1] * v2[0]
                                                            cp2 = v2[0] * v3[1] - v2[1] * v3[0]
                                                            if (cp1 * cp2 > 0):
                                                                newEl = np.array([n1,n2,n3,n4])
                                                            else:
                                                                newEl = np.array([n1,n2,n4,n3])
                                                            self.quadElements = np.array([[self.quadElements],[newEl]])
        
        newElMerged = elMerged
        newNodeElim = nodeElim
        return newElMerged,newNodeElim
        
        
    def mergeNestedTriEls(self, elMerged, el2El, nodeElim):
        """Object data modified: self.triElements
        Parameters
        ----------
        elMerged
        el2El
        nodeElim

        Returns
        -------
        newElMerged
        newNodeElim
        newEl2El
        """ 
        numEls = self.triElements.shape[0]
        for i in range(numEls):
            if (elMerged[i] == 0):
                for p in range(3):
                    j = el2El(i,p)
                    if (j > 0):
                        if (j != i and elMerged[i] == 0 and elMerged[j] == 0):
                            for q in range(3):
                                k = el2El(j,q)
                                if (k > 0):
                                    if (k != i and k != j and elMerged[i] == 0 and elMerged[j] == 0 and elMerged[k] == 0):
                                        allNds = np.array([self.triElements[i,:],self.triElements[j,:],self.triElements[k,:]])
                                        allNds = np.sort(allNds)
                                        commonNds = []
                                        for m in range(8):
                                            if (allNds[m] == allNds(m + 1)):
                                                commonNds = np.array([commonNds,allNds[m]])
                                        if (len(commonNds) == 5):
                                            center = 0
                                            for m in range(4):
                                                if (commonNds[m] == commonNds(m + 1)):
                                                    center = commonNds[m]
                                            if (center != 0):
                                                elMerged[i] = 1
                                                elMerged[j] = 1
                                                elMerged[k] = 1
                                                nodeElim[center] = 1
                                                newEl = []
                                                for m in range(5):
                                                    if (commonNds[m] != center):
                                                        newEl = np.array([newEl,commonNds[m]])
                                                self.triElements = np.array([[self.triElements],[newEl]])
                                                elMerged = np.array([[elMerged],[0]])
                                                allNbrs = np.array([el2El[i,:],el2El[j,:],el2El[k,:]])
                                                e2e = []
                                                for m in range(9):
                                                    nb = allNbrs[m]
                                                    if (nb != i and nb != j and nb != k):
                                                        e2e = np.array([e2e,nb])
                                                ln = len(e2e)
                                                if (ln != 3):
                                                    e2e = np.array([e2e,np.zeros((1,3 - ln))])
                                                el2El = np.array([[el2El],[e2e]])
        
        newElMerged = elMerged
        newNodeElim = nodeElim
        newEl2El = el2El
        return newElMerged,newNodeElim,newEl2El
        
        
    def mergeTriEls(self, elType): 
        """Object data modified: self.quadElements, self.nodes,
        self.triElements

        Parameters
        ----------
        elType

        Returns
        -------
        self
        """ 
        elMerged = np.zeros((self.triElements.shape[0],1))
        nodeElim = np.zeros((self.nodes.shape[0],1))
        el2El = np.zeros((self.triElements.shape[0],3))
        for i in np.arange(1,self.edges.shape[0]+1):
            if (self.edges(i,3) > 0 and self.edges(i,4) > 0):
                e1 = self.edges(i,3)
                e2 = self.edges(i,4)
                inserted = 0
                j = 1
                while (j <= 3 and inserted == 0):
    
                    if (el2El(e1,j) == 0):
                        el2El[e1,j] = e2
                        inserted = 1
                    else:
                        if (el2El(e1,j) == e2):
                            inserted = 1
                    j = j+1
    
                inserted = 0
                j = 1
                while (j <= 3 and inserted == 0):
    
                    if (el2El(e2,j) == 0):
                        el2El[e2,j] = e1
                        inserted = 1
                    else:
                        if (el2El(e2,j) == e1):
                            inserted = 1
                    j = j+1
    
        
        self.quadElements = []
        elMerged,nodeElim,self = self.mergeNestedQuadEls(elMerged,el2El,nodeElim)
        elMerged,nodeElim,el2El,self = self.mergeNestedTriEls(elMerged,el2El,nodeElim)
        if 'quad' in elType:
            for i in np.arange(1,self.triElements.shape[0]+1):
                if (elMerged[i] == 0):
                    sN1,ang = self.sortNodesByAngle(self.triElements[i,:])
                    for k in range(3):
                        j = el2El(i,k)
                        if (j > 0):
                            if (i != j and elMerged[j] == 0 and elMerged[i] == 0):
                                sN2,ang = self.sortNodesByAngle(self.triElements[j,:])
                                if ((sN1[1] == sN2[1] and sN1[2] == sN2[2]) or (sN1[1] == sN2[2] and sN1[2] == sN2[1])):
                                    elConn = np.array([sN1[0:2],sN2[0],sN1[2]])
                                    self.quadElements = np.array([[self.quadElements],[elConn]])
                                    elMerged[i] = 1
                                    elMerged[j] = 1
            elMerged,self = self.mergeElsBetweenAngles(elMerged,el2El,0.25 * np.pi,0.75 * np.pi)
            elMerged,self = self.mergeElsBetweenAngles(elMerged,el2El,0.16666 * np.pi,0.83333 * np.pi)
        
        newNds = []
        newLabel = np.zeros((self.nodes.shape[0],1))
        lab = 0
        for i in np.arange(1,self.nodes.shape[0]+1):
            if (nodeElim[i] == 0):
                lab = lab + 1
                newLabel[i] = lab
                newNds = np.array([[newNds],[self.nodes[i,:]]])
        
        self.nodes = newNds
        if ('quad' in elType):
            for i in range(self.triElements.shape[0]):
                if (elMerged[i] == 0):
                    elConn = np.array([self.triElements[i,:],0])
                    self.quadElements = np.array([[self.quadElements],[elConn]])
            for i in range(self.quadElements.shape[0]):
                for j in range(4):
                    orig = self.quadElements(i,j)
                    if (orig != 0):
                        self.quadElements[i,j] = newLabel(orig)
        else:
            newTriEls = []
            for i in range(self.triElements.shape[0]):
                if (elMerged[i] == 0):
                    newTriEls = np.array([[newTriEls],[self.triElements[i,:]]])
            for i in range(self.quadElements.shape[0]):
                nds = self.quadElements[i,:]
                triels = np.array([[nds[0],nds[1],nds[3]],[nds[2],nds[3],nds[1]]])
                newTriEls = np.array([[newTriEls],[triels]])
            self.triElements = newTriEls
            for i in range(self.triElements.shape[0]):
                for j in range(3):
                    orig = self.triElements(i,j)
                    if (orig != 0):
                        self.triElements[i,j] = newLabel(orig)
        
        return self
        
        
    def prepareForMesh(self):
        """Object data modified: self.maxEdgeLen, self.avgProjLen, self.maxElSize,
        self.MinElSize, self.nodeGL, self.edgeGL, self.triElGL
        Parameters
        ----------

        Returns
        -------
        self
        """ 
        numNds = self.nodes.shape[0]
        numEdges = self.edges.shape[0]
        self.triElements = np.array([])
        self.quadElements = np.array([])
        minSpacing = 10 * np.amax(np.abs(self.nodes))
        maxSpacing = 0
        avgSpacing = 0
        for i in range(numEdges):
            n1 = self.edges[i,0].astype(int)
            n2 = self.edges[i,1].astype(int)
            vec = self.nodes[n2,:] - self.nodes[n1,:]
            mag = np.linalg.norm(vec)
            avgSpacing = avgSpacing + mag
            if (mag < minSpacing):
                minSpacing = mag
            if (mag > maxSpacing):
                maxSpacing = mag
        
        self.maxEdgeLen = maxSpacing
        self.avgProjLen = 0.866025403784439 * (avgSpacing / self.edges.shape[0])
        self.maxElSize = 0.9 * maxSpacing
        self.minElSize = 0.8 * minSpacing
        minX = np.amin(self.nodes[:,0]) - maxSpacing
        maxX = np.amax(self.nodes[:,0]) + maxSpacing
        minY = np.amin(self.nodes[:,1]) - maxSpacing
        maxY = np.amax(self.nodes[:,1]) + maxSpacing
        avgSpacing = 0.5 * (minSpacing + maxSpacing)
        self.nodeGL = spatialGridList2D(minX,maxX,minY,maxY,avgSpacing,avgSpacing)
        for i in range(numNds):
            self.nodeGL = self.nodeGL.addEntry(i,self.nodes[i,:])
        
        self.edgeGL = spatialGridList2D(minX,maxX,minY,maxY,avgSpacing,avgSpacing)
        for i in range(numEdges):
            n1 = self.edges[i,0].astype(int)
            n2 = self.edges[i,1].astype(int)
            midPt = 0.5 * (self.nodes[n1,:] + self.nodes[n2,:])
            self.edgeGL = self.edgeGL.addEntry(i,midPt)
        
        self.triElGL = spatialGridList2D(minX,maxX,minY,maxY,avgSpacing,avgSpacing)
        self.getBoundaryEdgeNormals(minX,maxX,minY,maxY,np.sqrt(0.05) * minSpacing)
        #             gridSize = minSpacing/4;
    #             minX = min(self.nodes(:,1)) - 10*gridSize;
    #             maxX = max(self.nodes(:,1)) + 10*gridSize;
    #             minY = min(self.nodes(:,2)) - 10*gridSize;
    #             maxY = max(self.nodes(:,2)) + 10*gridSize;
    #             self.grdStat = gridStatus2D(minX,maxX,minY,maxY,gridSize);
    #             self.grdStat = self.grdStat.getInteriorGrid(self.nodes,self.edges);
    #             self.edges = self.grdStat.getBoundEdgeNormals(self.edges,self.nodes);
        return self
        
        
    def ptInEl(self, pt, elConn, nodes, marginFactor): 
        """Object data modified: none
        Parameters
        ----------
        pt
        elConn
        nodes
        marginFactor

        Returns
        -------
        inEl
        """ 
        inEl = 1
        n1 = elConn[0]
        n2 = elConn[1]
        n3 = elConn[2]
        elNd = np.array([[nodes[n1,:]],[nodes[n2,:]],[nodes[n3,:]]])
        centroid = 0.33333333 * (elNd[0,:] + elNd[1,:] + elNd[2,:])
        cenMat = np.ones((3,2)) * np.array([[centroid[0],0],[0,centroid[1]]])
        elNd = elNd - cenMat
        elNd = marginFactor * elNd
        elNd = elNd + cenMat
        v1 = elNd[1,:] - elNd[0,:]
        v2 = elNd[2,:] - elNd[0,:]
        v3 = pt - elNd[0,:]
        cp1 = v1[0] * v2[1] - v1[1] * v2[0]
        cp2 = v1[0] * v3[1] - v1[1] * v3[0]
        if (cp1 * cp2 <= 0):
            inEl = 0
        
        if (inEl == 1):
            v1 = elNd[2,:] - elNd[1,:]
            v2 = elNd[0,:] - elNd[1,:]
            v3 = pt - elNd[1,:]
            cp1 = v1[0] * v2[1] - v1[1] * v2[0]
            cp2 = v1[0] * v3[1] - v1[1] * v3[0]
            if (cp1 * cp2 <= 0):
                inEl = 0
        
        if (inEl == 1):
            v1 = elNd[0,:] - elNd[2,:]
            v2 = elNd[1,:] - elNd[2,:]
            v3 = pt - elNd[2,:]
            cp1 = v1[0] * v2[1] - v1[1] * v2[0]
            cp2 = v1[0] * v3[1] - v1[1] * v3[0]
            if (cp1 * cp2 <= 0):
                inEl = 0
        
        return inEl
        
        
    def skewNodes(self):
        """Object data modified: self.nodes
        Parameters
        ----------

        Returns
        -------
        self
        """ 
        skewMat = np.array([[1,np.tan(np.pi / 12)],[np.tan(np.pi / 12),1]])
        self.nodes = self.nodes @ skewMat
        return self
        
        
    def sortIncompleteEdges(self): 
        """Object data modified: self.edges
        Parameters
        ----------

        Returns
        -------
        self
        """ 
        completeEdges = np.array([])
        completeLabels = np.array([])
        incompleteEdges = np.array([])
        incompleteLabels = np.array([])
        for i in range(self.edges.shape[0]):
            if (self.edges[i,3] == 0):
                incompleteEdges = np.vstack([incompleteEdges,self.edges[i,:]]) if incompleteEdges.size else self.edges[[i],:]
                incompleteLabels = np.concatenate([incompleteLabels,[i]])
            else:
                completeEdges = np.vstack([completeEdges,[self.edges[i,:]]]) if completeEdges.size else self.edges[[i],:]
                completeLabels = np.concatenate([completeLabels,[i]])
        
        numInc = incompleteEdges.shape[0]
        eLen = np.array([])
        for i in range(numInc):
            n1 = incompleteEdges[i,0].astype(int)
            n2 = incompleteEdges[i,1].astype(int)
            vec = self.nodes[n1,:2] - self.nodes[n2,:2]
            mag = np.linalg.norm(vec)
            eLen = np.concatenate([eLen,[mag]])
        
        for i in range(numInc):
            for j in range(numInc-1):
                if (eLen[j+1] < eLen[j]):
                    edgeSwap = incompleteEdges[j,:]
                    incompleteEdges[j,:] = incompleteEdges[j+1,:]
                    incompleteEdges[j+1,:] = edgeSwap
                    swap = eLen[j]
                    eLen[j] = eLen[j+1]
                    eLen[j+1] = swap
                    labSwap = incompleteLabels[j]
                    incompleteLabels[j] = incompleteLabels[j+1]
                    incompleteLabels[j+1] = labSwap
        if incompleteEdges.size and completeEdges.size:
            self.edges = np.vstack([incompleteEdges,completeEdges])
            sortedLabels = np.concatenate([incompleteLabels,completeLabels]).astype(int)
        elif not incompleteEdges.size:
            self.edges = completeEdges
            sortedLabels = completeLabels.astype(int)
        elif not completeEdges.size:
            self.edges = incompleteEdges
            sortedLabels = incompleteLabels.astype(int)

        self.edgeGL = self.edgeGL.reOrderLabels(sortedLabels)
        return self
        
        
    def sortNodesByAngle(self, labels):
        """Object data modified: none
        Parameters
        ----------
        labels

        Returns
        -------
        sortedNodes
        angles
        """ 
        l1 = labels[0]
        l2 = labels[1]
        l3 = labels[2]
        angles = np.array([0,0,0])
        v1 = self.nodes[l2,:] - self.nodes[l1,:]
        mag1 = np.sqrt(v1 * v1.T)
        n1 = (1 / mag1) * v1
        v2 = self.nodes[l3,:] - self.nodes[l1,:]
        mag2 = np.sqrt(v2 * v2.T)
        n2 = (1 / mag2) * v2
        dp = n1 * n2.T
        angles[1] = np.arccos(dp)
        v1 = self.nodes[l1,:] - self.nodes[l2,:]
        mag1 = np.sqrt(v1 * v1.T)
        n1 = (1 / mag1) * v1
        v2 = self.nodes[l3,:] - self.nodes[l2,:]
        mag2 = np.sqrt(v2 * v2.T)
        n2 = (1 / mag2) * v2
        dp = n1 * n2.T
        angles[2] = np.arccos(dp)
        v1 = self.nodes[l1,:] - self.nodes[l3,:]
        mag1 = np.sqrt(v1 * v1.T)
        n1 = (1 / mag1) * v1
        v2 = self.nodes[l2,:] - self.nodes[l3,:]
        mag2 = np.sqrt(v2 * v2.T)
        n2 = (1 / mag2) * v2
        dp = n1 * n2.T
        angles[3] = np.arccos(dp)
        sortedNodes = labels
        for i in range(3):
            for j in range(2):
                if (angles[j+1] > angles[j]):
                    swap = sortedNodes[j]
                    sortedNodes[j] = sortedNodes[j+1]
                    sortedNodes[j+1] = swap
                    swap = angles[j]
                    angles[j] = angles[j+1]
                    angles[j+1] = swap
        
        return sortedNodes,angles
        
        
    def uniformBoundarySpacing(self,maxIt):
        """Object data modified: self.nodes

        Parameters
        ----------
        maxIt : int

        Returns
        ---------
        self
        
        """
        numNds = self.nodes.shape[0]
        numEdges = self.edges.shape[0]
        eqnFact = 10
        cols = 2 * numNds
        rows = 2 * cols
        eqnMat = np.zeros((rows,cols))
        for i in range(numNds):
            j = (2 * i) + 1
            eqnMat[j-1,j-1] = 1
            eqnMat[j,j] = 1
        
        for i in range(numEdges):
            n1 = self.edges[i,0].astype(int)
            n2 = self.edges[i,1].astype(int)
            j = cols + (2 * i) + 1
            k11 = (2 * n1)
            k12 = (2 * n1) + 1
            k21 = (2 * n2)
            k22 = (2 * n2) + 1
            eqnMat[j-1,k21] = eqnFact
            eqnMat[j-1,k11] = - eqnFact
            eqnMat[j,k22] = eqnFact
            eqnMat[j,k12] = - eqnFact
        
        Q,R = np.linalg.qr(eqnMat)
        for it in range(maxIt):
            edgeDir = self.getBoundaryData()
            rhs = np.zeros((rows,1))
            for i in range(numNds):
                j = (2 * i) + 1
                rhs[j-1] = self.nodes[i,0]
                rhs[j] = self.nodes[i,1]
            for i in range(numEdges):
                j = cols + (2 * i) + 1
                rhs[j-1] = eqnFact * edgeDir[i,0]
                rhs[j] = eqnFact * edgeDir[i,1]
            solnVec = np.linalg.solve(R,Q.T) @ rhs
            self.nodes = np.array([])
            for i in range(numNds):
                j = (2 * i) + 1
                nd = solnVec[j-1:j+1].T
                self.nodes = np.vstack([self.nodes,nd]) if self.nodes.size else nd
        
        return self
        
        
    def unSkewNodes(self = None): 
        skewMat = np.array([[1,np.tan(np.pi / 12)],[np.tan(np.pi / 12),1]])
        invSkew = np.linalg.solve(skewMat,np.array([[1,0],[0,1]]))
        self.nodes = self.nodes * invSkew
        return self
        

class NuMesh3D():
    #UNTITLED Summary of this class goes here
    #   Detailed explanation goes here
    
    def __init__(self,boundaryNodes = None,boundaryFaces = None): 
        self.nodes = []
        self.elements = []
        self.boundaryFaces = []
        self.nodes = boundaryNodes
        self.boundaryFaces = boundaryFaces
        
        
    def createSweptMesh(self,sweepMethod = None,direction = None,sweepDistance = None,sweepElements = None,followNormal = None,destNodes = None): 
        if ('to_dest_nodes' in sweepMethod):
            totalSwpEls = np.sum(sweepElements, 'all')
            guideNodes = []
            for j in range(self.nodes.shape[0]):
                guideNodes = np.array([[guideNodes],[np.transpose(self.nodes[j,:])]])
            destLen = destNodes.shape[1]
            for i in range(destLen):
                Dnds = destNodes[i]
                col = []
                for j in range(self.nodes.shape[0]):
                    col = np.array([[col],[np.transpose(Dnds[j,:])]])
                guideNodes = np.array([guideNodes,col])
            allNodes = []
            guideParVal = 0
            for i in range(len(sweepElements)):
                nextEnt = guideParVal[-1] + sweepElements[i]
                guideParVal = np.array([guideParVal,nextEnt])
            guideParVal = (1 / nextEnt) * guideParVal
            #                 guideParVal = linspace(0,1,destLen+1);
            allParVal = np.linspace(0,1,totalSwpEls + 1)
            for i in range(guideNodes.shape[0]):
                row = interpolator_wrap(guideParVal,guideNodes[i,:],allParVal,'pchip')
                allNodes = np.array([[allNodes],[row]])
            nodes = []
            for i in range(totalSwpEls+1):
                for j in range(0,allNodes.shape[0],3):
                    nodes = np.array([[nodes],[np.transpose(allNodes[j:j+2,i])]])
            numLayerNds = self.nodes.shape[0]
            numLayerEls = self.boundaryFaces.shape[0]
            elements = []
            four1s = np.ones((1,4))
            for i in range(totalSwpEls):
                for j in range(numLayerEls):
                    face = self.boundaryFaces[j,:]
                    lowFace = face + numLayerNds * (i-1) * four1s
                    highFace = face + numLayerNds * i * four1s
                    if (face(4) == 0):
                        newEl = np.array([lowFace[0,2],highFace[0,2],0,0])
                        n1 = newEl[0]
                        n2 = newEl[1]
                        n3 = newEl[2]
                        n4 = newEl[3]
                        v1 = nodes[n2,:] - nodes[n1,:]
                        v2 = nodes[n3,:] - nodes[n1,:]
                        v3 = nodes[n4,:] - nodes[n1,:]
                        vMat = np.array([[v1],[v2],[v3]])
                        detV = np.linalg.det(vMat)
                        if (detV < 0):
                            swap = newEl[1]
                            newEl[1] = newEl[2]
                            newEl[2] = swap
                            swap = newEl[4]
                            newEl[4] = newEl[5]
                            newEl[5] = swap
                    else:
                        newEl = np.array([lowFace,highFace])
                        n1 = newEl[0]
                        n2 = newEl[1]
                        n3 = newEl[2]
                        n5 = newEl[4]
                        v1 = nodes[n2,:] - nodes[n1,:]
                        v2 = nodes[n3,:] - nodes[n1,:]
                        v3 = nodes[n5,:] - nodes[n1,:]
                        vMat = np.array([[v1],[v2],[v3]])
                        detV = np.linalg.det(vMat)
                        if (detV < 0):
                            swap = newEl[1]
                            newEl[1] = newEl[3]
                            newEl[3] = swap
                            swap = newEl[5]
                            newEl[5] = newEl[7]
                            newEl[7] = swap
                    elements = np.array([[elements],[newEl]])
            self.nodes = nodes
            self.elements = elements
        
        return nodes,elements,self
        

class shellRegion: 
    """
    Attributes
    -----------
    type : str
    keyPts : list
    edgeEls : list
    """
    def __init__(self, regType, keyPoints, numEdgeEls): 

        self.type = regType
        self.keyPts = keyPoints
        self.edgeEls = numEdgeEls
        
    
    def createShellMesh(self, elType, method):
        """Object data modified: none
        Parameters
        ----------
        elType
        method : str
        
        Returns
        -------
        nodes
        elements
        """
        if 'structured' in method:
            if ('quad' in self.type) or ('cyl' in self.type):
                xNodes = np.amax(self.edgeEls[[0,2]]) + 1
                yNodes = np.amax(self.edgeEls[[1,3]]) + 1
                boundaryNodes = np.stack([np.linspace(- 1,1,xNodes),- 1 * np.ones((xNodes))],axis=1)
                boundaryEdges = np.stack([range(xNodes-1),range(1,xNodes)],axis=1)
                mesh = NuMesh2D(boundaryNodes,boundaryEdges)
                nodes,elements = mesh.createSweptMesh('in_direction',np.array([0,1]),2,(yNodes - 1),0,[])
                ndElim = np.zeros((len(nodes),1)).astype(int)

                if (self.edgeEls[0] < self.edgeEls[2]):
                    nds2Elim = (self.edgeEls[2] - self.edgeEls[0])
                    rowstep = np.ceil((xNodes - nds2Elim) / (nds2Elim + 1)).astype(int) + 1
                    for i in range(rowstep,xNodes+rowstep,rowstep):
                        nd = i
                        ndElim[nd] = nd - 1
                    xInc = 2 / self.edgeEls[0]
                    xPrev = - 1
                    for i in range(1,xNodes):
                        nd = i
                        if (ndElim[nd] == 0):
                            nodes[nd,0] = xPrev + xInc
                            xPrev = xPrev + xInc

                elif (self.edgeEls[2] < self.edgeEls[0]):
                    nds2Elim = (self.edgeEls[0] - self.edgeEls[2])
                    rowstep = np.ceil((xNodes - nds2Elim) / (nds2Elim + 1)).astype(int) + 1
                    for i in range(rowstep,xNodes,rowstep):
                        nd = (yNodes - 1) * xNodes + i
                        ndElim[nd] = nd - 1
                    xInc = 2 / self.edgeEls[2]
                    xPrev = - 1
                    for i in range(1,xNodes):
                        nd = (yNodes - 1) * xNodes + i
                        if (ndElim[nd] == 0):
                            nodes[nd,0] = xPrev + xInc
                            xPrev = xPrev + xInc

                if (self.edgeEls[1] < self.edgeEls[3]):
                    nds2Elim = (self.edgeEls[3] - self.edgeEls[1])[0]
                    rowstep = np.ceil((yNodes - nds2Elim) / (nds2Elim + 1)).astype(int) + 1
                    for j in range(rowstep,yNodes+rowstep,rowstep):
                        nd = j * xNodes
                        ndElim[nd] = nd - xNodes
                    yInc = 2 / self.edgeEls[1]
                    yPrev = - 1
                    for j in range(1,yNodes):
                        nd = j * xNodes
                        if (ndElim[nd] == 0):
                            nodes[nd,1] = yPrev + yInc
                            yPrev = yPrev + yInc

                elif (self.edgeEls[3] < self.edgeEls[1]):
                    nds2Elim = (self.edgeEls[1] - self.edgeEls[3])[0]
                    rowstep = np.ceil((yNodes - nds2Elim) / (nds2Elim + 1)).astype(int) + 1
                    for j in range(rowstep,yNodes+rowstep,rowstep):
                        nd = (j - 1) * xNodes + 1
                        ndElim[nd] = nd - xNodes
                    yInc = 2 / self.edgeEls[3]
                    yPrev = - 1
                    for j in range(1,yNodes):
                        nd = (j - 1) * xNodes + 1
                        if (ndElim[nd] == 0):
                            nodes[nd,1] = yPrev + yInc
                            yPrev = yPrev + yInc

                if (np.sum(ndElim) > 0.1):
                    newNds = np.ndarray([])
                    newLabel = np.zeros((len(nodes)))
                    lab = 0
                    for i in range(len(nodes)):
                        if (ndElim[i] == 0):
                            lab = lab + 1
                            newLabel[i] = lab
                            try:
                                newNds = np.vstack([newNds,nodes[i,:]])
                            except ValueError:
                                newNds = nodes[i,:]
                    nodes = newNds
                    for i in range(len(elements)):
                        newEl = np.array([])
                        for j in range(4):
                            orig = elements[i,j]
                            if (orig != 0):
                                if (ndElim[orig] != 0):
                                    nlab = newLabel[ndElim[orig]]
                                else:
                                    nlab = newLabel[orig]
                                if (nlab != 0):
                                    if (not np.any(newEl == nlab) ):
                                        try:
                                            newEl = np.hstack([newEl,nlab])
                                        except ValueError:
                                            newEl = nlab
                        ln = len(newEl)
                        if (ln < 4):
                            newEl = np.hstack([newEl,np.zeros((4 - ln))])
                        elements[i,:] = newEl.reshape(-1)
                newNodes = np.array([])

                for i in range(len(nodes)):
                    XYZ = self.XYZCoord(nodes[i,:])
                    try:
                        newNodes = np.vstack([newNodes,XYZ])
                    except ValueError:
                        newNodes = XYZ
                nodes = newNodes

            else:
                raise Exception('Only quadrilateral or cylinder shell regions can use the structured meshing option')

        else:
            boundaryNodes,edges = self.initialBoundary()
            # plot2DMesh(boundaryNodes,[])
            # keyboard
            mesh = NuMesh2D(boundaryNodes,edges[:,:2])
            etaNodes,elements,mesh = mesh.createPlanarMesh(elType,1)
            nodes = []
            for i in range(len(etaNodes)):
                XYZ = self.XYZCoord(etaNodes[i,:])
                nodes = np.array([[nodes],[XYZ]])
        
        return nodes,elements
        
        
    def initialBoundary(self):
        """ Object data modified: none
        Parameters
        ----------

        Returns
        -------
        nodes
        edges
        """
        if 'quad' in self.type:
            nodes = np.array([])
            edges = np.array([])
            ## initialize nodes on side 1
            delE = 2 / self.edgeEls[0]
            for i in range(1, self.edgeEls[0]+1):
                e1 = delE * i - 1
                pt = np.array([e1,-1])
                nodes = np.vstack([nodes,pt]) if nodes.size else pt
            ## initialize nodes on side 2
            delE = 2 / self.edgeEls[1]
            for i in range(1, self.edgeEls[1]+1):
                e2 = delE * i - 1
                pt = np.array([1,e2])
                nodes = np.vstack([nodes,pt])
            ## initialize nodes on side 3
            delE = 2 / self.edgeEls[2]
            for i in range(1, self.edgeEls[2]+1):
                e1 = 1 - delE * i
                pt = np.array([e1,1])
                nodes = np.vstack([nodes,pt])
            ## initialize nodes on side 4
            delE = 2 / self.edgeEls[3]
            for i in range(1, self.edgeEls[3]+1):
                e2 = 1 - delE * i
                pt = np.array([-1,e2])
                nodes = np.vstack([nodes,pt])
            ## create edges
            numNds = len(nodes)
            edges = np.concatenate([
                np.array([range(numNds)]).T,
                np.concatenate([np.array([range(1,numNds)]).T, [[0]]],axis=0),
                - np.ones((numNds,1)),np.zeros((numNds,3))
            ],axis=1).astype(int)
        elif 'tri' in self.type:
            nodes = []
            ## initialize nodes on side 1
            delE = 1 / self.edgeEls[0]
            for i in range(1, self.edgeEls[0]+1):
                e1 = delE * i
                pt = np.array([e1,0])
                nodes = np.concatenate([nodes,pt])
            ## initialize nodes on side 2
            delE = 1 / self.edgeEls[1]
            for i in range(1, self.edgeEls[1]+1):
                e1 = 1 - delE * i
                e2 = delE * i
                pt = np.array([e1,e2])
                nodes = np.concatenate([nodes,pt])
            ## initialize nodes on side 3
            delE = 1 / self.edgeEls[2]
            for i in range(1, self.edgeEls[2]+1):
                e2 = 1 - delE * i
                pt = np.array([0,e2])
                nodes = np.concatenate([nodes,pt])
            ## create edges
            numNds = len(nodes)
            edges = np.array([np.transpose(np.array([range(numNds+1)])),np.transpose(np.array([range(1,numNds),1])),- np.ones((numNds,1)),np.zeros((numNds,3))])
        elif 'sphere' in self.type:
            nodes = []
            delTheta = 2 * np.pi / self.edgeEls[0]
            for i in range(1, self.edgeEls[0]+1):
                theta = delTheta * i
                x = 0.5 * np.pi * np.cos(theta)
                y = 0.5 * np.pi * np.sin(theta)
                nodes = np.concatenate([nodes,np.array([x,y])])
            numNds = len(nodes)
            edges = np.array([np.transpose(np.array([range(numNds)])),np.transpose(np.array([range(1,numNds),1])),- np.ones((numNds,1)),np.zeros((numNds,3))])

        return nodes,edges
        
        
    def XYZCoord(self, eta): 
        """
        Parameters
        ----------
        eta

        Returns
        -------
        XYZ
        """
        if 'quad4' == self.type:
            Nvec = np.zeros((1,4))
            Nvec[0] = 0.25 * (eta[0] - 1) * (eta[1] - 1)
            Nvec[1] = - 0.25 * (eta[0] + 1) * (eta[1] - 1)
            Nvec[2] = 0.25 * (eta[0] + 1) * (eta[1] + 1)
            Nvec[3] = - 0.25 * (eta[0] - 1) * (eta[1] + 1)
            XYZ = Nvec * self.keyPts
        elif 'quad9' == self.type:
            r1 = - 1
            r2 = 0
            r3 = 1
            Nvec = np.zeros((1,9))
            Nvec[0] = 0.25 * (eta[0] - r2) * (eta[0] - r3) * (eta[1] - r2) * (eta[1] - r3)
            Nvec[1] = 0.25 * (eta[0] - r1) * (eta[0] - r2) * (eta[1] - r2) * (eta[1] - r3)
            Nvec[2] = 0.25 * (eta[0] - r1) * (eta[0] - r2) * (eta[1] - r1) * (eta[1] - r2)
            Nvec[3] = 0.25 * (eta[0] - r2) * (eta[0] - r3) * (eta[1] - r1) * (eta[1] - r2)
            Nvec[4] = - 0.5 * (eta[0] - r1) * (eta[0] - r3) * (eta[1] - r2) * (eta[1] - r3)
            Nvec[5] = - 0.5 * (eta[0] - r1) * (eta[0] - r2) * (eta[1] - r1) * (eta[1] - r3)
            Nvec[6] = - 0.5 * (eta[0] - r1) * (eta[0] - r3) * (eta[1] - r1) * (eta[1] - r2)
            Nvec[7] = - 0.5 * (eta[0] - r2) * (eta[0] - r3) * (eta[1] - r1) * (eta[1] - r3)
            Nvec[8] = (eta[0] - r1) * (eta[0] - r3) * (eta[1] - r1) * (eta[1] - r3)
            XYZ = Nvec * self.keyPts
        elif 'quad16' == self.type:
            r1 = - 1
            r2 = - 0.333333333333333
            r3 = 0.333333333333333
            r4 = 1
            coef = np.array([0.31640625,- 0.31640625,0.31640625,- 0.31640625,- 0.94921875,0.94921875,0.94921875,- 0.94921875,- 0.94921875,0.94921875,0.94921875,- 0.94921875,2.84765625,- 2.84765625,2.84765625,- 2.84765625])
            Nvec = np.zeros((16))
            Nvec[0] = (eta[0] - r2) * (eta[0] - r3) * (eta[0] - r4) * (eta[1] - r2) * (eta[1] - r3) * (eta[1] - r4)
            Nvec[1] = (eta[0] - r1) * (eta[0] - r2) * (eta[0] - r3) * (eta[1] - r2) * (eta[1] - r3) * (eta[1] - r4)
            Nvec[2] = (eta[0] - r1) * (eta[0] - r2) * (eta[0] - r3) * (eta[1] - r1) * (eta[1] - r2) * (eta[1] - r3)
            Nvec[3] = (eta[0] - r2) * (eta[0] - r3) * (eta[0] - r4) * (eta[1] - r1) * (eta[1] - r2) * (eta[1] - r3)
            Nvec[4] = (eta[0] - r1) * (eta[0] - r3) * (eta[0] - r4) * (eta[1] - r2) * (eta[1] - r3) * (eta[1] - r4)
            Nvec[5] = (eta[0] - r1) * (eta[0] - r2) * (eta[0] - r4) * (eta[1] - r2) * (eta[1] - r3) * (eta[1] - r4)
            Nvec[6] = (eta[0] - r1) * (eta[0] - r2) * (eta[0] - r3) * (eta[1] - r1) * (eta[1] - r3) * (eta[1] - r4)
            Nvec[7] = (eta[0] - r1) * (eta[0] - r2) * (eta[0] - r3) * (eta[1] - r1) * (eta[1] - r2) * (eta[1] - r4)
            Nvec[8] = (eta[0] - r1) * (eta[0] - r2) * (eta[0] - r4) * (eta[1] - r1) * (eta[1] - r2) * (eta[1] - r3)
            Nvec[9] = (eta[0] - r1) * (eta[0] - r3) * (eta[0] - r4) * (eta[1] - r1) * (eta[1] - r2) * (eta[1] - r3)
            Nvec[10] = (eta[0] - r2) * (eta[0] - r3) * (eta[0] - r4) * (eta[1] - r1) * (eta[1] - r2) * (eta[1] - r4)
            Nvec[11] = (eta[0] - r2) * (eta[0] - r3) * (eta[0] - r4) * (eta[1] - r1) * (eta[1] - r3) * (eta[1] - r4)
            Nvec[12] = (eta[0] - r1) * (eta[0] - r3) * (eta[0] - r4) * (eta[1] - r1) * (eta[1] - r3) * (eta[1] - r4)
            Nvec[13] = (eta[0] - r1) * (eta[0] - r2) * (eta[0] - r4) * (eta[1] - r1) * (eta[1] - r3) * (eta[1] - r4)
            Nvec[14] = (eta[0] - r1) * (eta[0] - r2) * (eta[0] - r4) * (eta[1] - r1) * (eta[1] - r2) * (eta[1] - r4)
            Nvec[15] = (eta[0] - r1) * (eta[0] - r3) * (eta[0] - r4) * (eta[1] - r1) * (eta[1] - r2) * (eta[1] - r4)
            Nvec = np.multiply(coef,Nvec)
            XYZ = np.matmul(Nvec,self.keyPts)
        elif 'tri3' == self.type:
            Nvec = np.zeros((1,3))
            Nvec[0] = 1 - eta[0] - eta[1]
            Nvec[1] = eta[0]
            Nvec[2] = eta[1]
            XYZ = Nvec * self.keyPts
        elif 'tri6' == self.type:
            Nvec = np.zeros((1,6))
            Nvec[0] = 2 * (eta[0] + eta[1] - 1) * (eta[0] + eta[1] - 0.5)
            Nvec[1] = 2 * eta[0] * (eta[0] - 0.5)
            Nvec[2] = 2 * eta[1] * (eta[1] - 0.5)
            Nvec[3] = - 4 * eta[0] * (eta[0] + eta[1] - 1)
            Nvec[4] = 4 * eta[0] * eta[1]
            Nvec[5] = - 4 * eta[1] * (eta[0] + eta[1] - 1)
            XYZ = np.matmul(Nvec,self.keyPts)
        elif 'tri10' == self.type:
            r2 = 1 / 3
            r3 = 2 / 3
            coef = np.array([- 4.5,4.5,4.5,13.5,- 13.5,13.5,13.5,- 13.5,13.5,- 27])
            Nvec = np.zeros((1,10))
            Nvec[0] = (eta[0] + eta[1] - r2) * (eta[0] + eta[1] - r3) * (eta[0] + eta[1] - 1)
            Nvec[1] = eta[0] * (eta[0] - r2) * (eta[0] - r3)
            Nvec[2] = eta[1] * (eta[1] - r2) * (eta[1] - r3)
            Nvec[3] = eta[0] * (eta[0] + eta[1] - r3) * (eta[0] + eta[1] - 1)
            Nvec[4] = eta[0] * (eta[0] - r2) * (eta[0] + eta[1] - 1)
            Nvec[5] = eta[0] * eta[1] * (eta[0] - r2)
            Nvec[6] = eta[0] * eta[1] * (eta[1] - r2)
            Nvec[7] = eta[1] * (eta[1] - r2) * (eta[0] + eta[1] - 1)
            Nvec[8] = eta[1] * (eta[0] + eta[1] - r3) * (eta[0] + eta[1] - 1)
            Nvec[9] = eta[0] * eta[1] * (eta[0] + eta[1] - 1)
            Nvec = np.multiply(coef,Nvec)
            XYZ = np.matmul(Nvec,self.keyPts)

        elif 'sphere' == self.type:
            vec = self.keyPts[1,:] - self.keyPts[0,:]
            outerRad = np.linalg.norm(vec)
            phiComp = np.sqrt(eta * eta.T)
            if (eta[0] > 0):
                theta = np.arctan(eta[1] / eta[0])
            else:
                if (eta[0] < 0):
                    theta = np.pi + np.arctan(eta[1] / eta[0])
                else:
                    theta = np.arctan(eta[1] / 1e-08)
            phi = 0.5 * np.pi - phiComp
            xloc = outerRad * np.cos(theta) * np.cos(phi)
            yloc = outerRad * np.sin(theta) * np.cos(phi)
            zloc = outerRad * np.sin(phi)
            a1 = (1 / outerRad) * vec
            vec2 = self.keyPts[2,:] - self.keyPts[1,:]
            vec3 = np.array([
                (vec[1] * vec2[2] - vec[2] * vec2[1]),
                (vec[2] * vec2[0] - vec[0] * vec2[2]),
                (vec[0] * vec2[1] - vec[1] * vec2[0])
                ])
            mag = np.sqrt(vec3 * vec3.T)
            a3 = (1 / mag) * vec3
            a2 = np.array([
                (a3[1] * a1[2] - a3[2] * a1[1]),
                (a3[2] * a1[0] - a3[0] * a1[2]),
                (a3[0] * a1[1] - a3[1] * a1[0])
                ])
            alpha = np.array([[a1],[a2],[a3]])
            XYZ = (np.array([xloc,yloc,zloc]) * alpha + self.keyPts[1,:])

        return XYZ


class elementSet:
    """
    Attributes
    ----------
    name : str
    plygroups : list[Plygroup]
    elementsList : array
    """
    def __init__(self, setName, setPlyGroups, setElList):
        self.name = setName
        self.plygroups = setPlyGroups
        self.elementList = setElList


class spatialGridList2D():
    """

    Attributes
    ----------
    firstEnt : array
    label : array
    nextEnt : array
    xMin : int
    yMin : int
    xGSz : int
        x grid size
    yGSz : int
        y grid size
    xRows : int
    yRows : int
    """

    def __init__(
            self,
            minimumX, 
            maximumX,
            minimumY,
            maximumY,
            xGridSize,
            yGridSize
        ): 
        self.firstEnt = np.array([], dtype=int)
        self.label = np.array([], dtype=int)
        self.nextEnt = np.array([], dtype=int)
        self.xMin = 0
        self.yMin = 0
        self.xGSz = 1
        self.yGSz = 1
        self.xRows = 1
        self.yRows = 1
        self.xMin = minimumX
        self.yMin = minimumY
        self.xGSz = xGridSize
        self.yGSz = yGridSize
        XLen = maximumX - minimumX
        YLen = maximumY - minimumY
        self.xRows = np.ceil(XLen / xGridSize).astype(int)
        self.yRows = np.ceil(YLen / yGridSize).astype(int)
        # NOTE: changed firstEnt from array of 0s to array of -1
        # since 0 is used as a valid index in python and
        # should be distinguished from default value.
        self.firstEnt = np.full((self.xRows,self.yRows),-1).astype(int)
        return
        
        
    def addEntry(self, val, coord): 
        """
        Object variables modified: self.label, self.nextEnt, self.firstEnt

        Parameters
        ----------
        val
        coord

        Returns
        -------
        self
        
        """
        xRow = np.ceil((coord[0] - self.xMin) / self.xGSz).astype(int) - 1
        yRow = np.ceil((coord[1] - self.yMin) / self.xGSz).astype(int) - 1
        if (self.firstEnt[xRow,yRow] == -1):
            self.label = np.concatenate([self.label,[val]])
            self.nextEnt = np.concatenate([self.nextEnt,[0]])
            self.firstEnt[xRow,yRow] = len(self.label) - 1
        else:
            i = self.firstEnt[xRow,yRow]
            inserted = 0
            while (inserted == 0):
    
                if (self.nextEnt[i] == 0):
                    self.label = np.concatenate([self.label,[val]])
                    self.nextEnt = np.concatenate([self.nextEnt,[0]])
                    self.nextEnt[i] = len(self.label) - 1
                    inserted = 1
                else:
                    i = self.nextEnt[i]
    
        return self
        
        
    def findInRadius(self, point, radius): 
        """
        Parameters
        ----------
        point
        radius
        """
        labelList = self.findInXYMargin(point, radius, radius)
        return labelList
        
        
    def findInXYMargin(self, point, Xmargin, Ymargin): 
        """Object variables modified: none
        Parameters
        ----------
        point
        Xmargin
        Ymargin
        """
        if (Xmargin == -1):
            iMax = self.xRows
            iMin = 1
        else:
            iMax = np.ceil((point[0] + Xmargin - self.xMin) / self.xGSz).astype(int)
            iMax = np.amin(np.array([iMax,self.xRows]))
            iMin = np.ceil((point[0] - Xmargin - self.xMin) / self.xGSz).astype(int)
            iMin = np.amax(np.array([iMin,1]))
        
        if (Ymargin == - 1):
            jMax = self.yRows
            jMin = 1
        else:
            jMax = np.ceil((point[1] + Ymargin - self.yMin) / self.yGSz).astype(int)
            jMax = np.amin(np.array([jMax,self.yRows]))
            jMin = np.ceil((point[1] - Ymargin - self.yMin) / self.yGSz).astype(int)
            jMin = np.amax(np.array([jMin,1]))
        
        # adjust for python indexing
        iMin += -1
        iMax += -1
        jMin += -1
        jMax += -1

        labelList = np.array([],dtype=int)
        for i in range(iMin,iMax+1):
            for j in range(jMin,jMax+1):
                if (self.firstEnt[i,j] != -1):
                    k = self.firstEnt[i,j]
                    endReached = 0
                    while (endReached == 0):  
                        labelList = np.concatenate([labelList,[self.label[k]]])
                        if (self.nextEnt[k] == 0):
                            endReached = 1
                        else:
                            k = self.nextEnt[k]
    
        return labelList
        
        
    def reOrderLabels(self, labelOrder):
        """Object variables modified: self.label
        Parameters
        ----------
        labelOrder
        """
        currentLabel = np.zeros((np.amax(labelOrder)+1))
        for i in range(len(labelOrder)):
            currentLabel[labelOrder[i]] = i
        
        for i in range(len(self.label)):
            orig = self.label[i]
            if (orig != 0):
                self.label[i] = currentLabel[orig]
        
        return self
        

class spatialGridList3D():
    """
    Parameters
    ----------
    firstEnt
    label
    nextEnt
    xMin
    yMin
    zMin
    xGSz
    yGSz
    zGSz
    xRow
    yRow
    zRow
    """
    

        
    def __init__(self,minimumX,maximumX,minimumY,maximumY,minimumZ,maximumZ,xGridSize,yGridSize,zGridSize): 
        self.firstEnt = np.array([])
        self.label = np.array([], dtype=int)
        self.nextEnt = np.array([])
        self.xMin = 0
        self.yMin = 0
        self.zMin = 0
        self.xGSz = 1
        self.yGSz = 1
        self.zGSz = 1
        self.xRows = 1
        self.yRows = 1
        self.zRows = 1
        self.xMin = minimumX
        self.yMin = minimumY
        self.zMin = minimumZ
        self.xGSz = xGridSize
        self.yGSz = yGridSize
        self.zGSz = zGridSize
        XLen = maximumX - minimumX
        YLen = maximumY - minimumY
        ZLen = maximumZ - minimumZ
        self.xRows = np.ceil(XLen / xGridSize).astype(int)
        self.yRows = np.ceil(YLen / yGridSize).astype(int)
        self.zRows = np.ceil(ZLen / zGridSize).astype(int)
        self.firstEnt = np.zeros((self.xRows,self.yRows,self.zRows)).astype(int)
        return
        
        
    def addEntry(self, val, coord): 
        """
        """
        xRow = np.ceil((coord[0] - self.xMin) / self.xGSz).astype(int) - 1
        yRow = np.ceil((coord[1] - self.yMin) / self.yGSz).astype(int) - 1
        zRow = np.ceil((coord[2] - self.zMin) / self.zGSz).astype(int) - 1
        if (self.firstEnt[xRow,yRow,zRow] == 0):
            self.label = np.concatenate([self.label,[val]]).astype(int)
            self.nextEnt = np.concatenate([self.nextEnt,[0]]).astype(int)
            self.firstEnt[xRow,yRow,zRow] = len(self.label) - 1
        else:
            i = self.firstEnt[xRow,yRow,zRow].astype(int)
            inserted = 0
            while (inserted == 0):
    
                if (self.nextEnt[i] == 0):
                    self.label = np.concatenate([self.label,[val]])
                    self.nextEnt = np.concatenate([self.nextEnt,[0]])
                    self.nextEnt[i] = len(self.label) - 1
                    inserted = 1
                else:
                    i = self.nextEnt[i]
    
        
        return self
        
        
    def findInRadius(self,point,radius): 
        labelList = self.findInXYZMargin(point,radius,radius,radius)
        return labelList
        
        
    def findInXYZMargin(self,point,Xmargin,Ymargin,Zmargin): 
        """ Object data modified: none

        Parameters
        ----------
        point
        Xmargin
        Ymargin
        Zmargin

        Returns
        -------
        labelList
        """
        if (Xmargin == - 1):
            iMax = self.xRows
            iMin = 1
        else:
            iMax = np.ceil((point[0] + Xmargin - self.xMin) / self.xGSz).astype(int)
            iMax = np.amin(np.array([iMax,self.xRows]))
            iMin = np.ceil((point[0] - Xmargin - self.xMin) / self.xGSz).astype(int)
            iMin = np.amax(np.array([iMin,1]))
        
        if (Ymargin == - 1):
            jMax = self.yRows
            jMin = 1
        else:
            jMax = np.ceil((point[1] + Ymargin - self.yMin) / self.yGSz).astype(int)
            jMax = np.amin(np.array([jMax,self.yRows]))
            jMin = np.ceil((point[1] - Ymargin - self.yMin) / self.yGSz).astype(int)
            jMin = np.amax(np.array([jMin,1]))
        
        if (Zmargin == - 1):
            kMax = self.zRows
            kMin = 1
        else:
            kMax = np.ceil((point[2] + Zmargin - self.zMin) / self.zGSz).astype(int)
            kMax = np.amin(np.array([kMax,self.zRows]))
            kMin = np.ceil((point[2] - Zmargin - self.zMin) / self.zGSz).astype(int)
            kMin = np.amax(np.array([kMin,1]))
        
        # adjust for python indexing
        iMin += -1
        iMax += -1
        jMin += -1
        jMax += -1
        kMin += -1
        kMax += -1

        labelList = np.array([],dtype=int)
        for i in range(iMin,iMax+1):
            for j in range(jMin,jMax+1):
                for k in range(kMin,kMax+1):
                    if (self.firstEnt[i,j,k] != 0):
                        m = self.firstEnt[i,j,k]
                        endReached = 0
                        while (endReached == 0):
                            labelList = np.concatenate([labelList,[self.label[m]]])
                            if (self.nextEnt[m] == 0):
                                endReached = 1
                            else:
                                m = self.nextEnt[m]
    
        return labelList
        
        
    def reOrderLabels(self,labelOrder): 
        """Object data modified: self.label

        Parameters
        ----------
        labelOrder

        Returns
        -------
        self
        """

        # I think something like this would work and is a little more concise:
        # currentLabel = [self.label[i] for i in labelOrder]
        # self.label = currentLabel
        currentLabel = np.zeros((np.amax(labelOrder),1))
        for i in range(len(labelOrder)):
            currentLabel[labelOrder[i]] = i
        
        for i in range(len(self.label)):
            orig = self.label[i]
            if (orig != 0):
                self.label[i] = currentLabel[orig]
        
        return self

class edges():
    def __init__(self, coords):
        self.coords = coords
        # self.otherstuff = other stuff