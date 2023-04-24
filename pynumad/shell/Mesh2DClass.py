import numpy as np
from pynumad.shell.SpatialGridList2DClass import *
import plotly.graph_objects as go

class Mesh2D():

    def __init__(self,boundaryNodes,boundaryEdges=[]):
        self.nodeGL = None
        self.edgeGL = None
        self.triElGL = None
        
        self.minEdgeLen = 0.0
        self.maxEdgeLen = 1.0
        self.avgProjLen = 1.0
        
        self.numBndNodes = len(boundaryNodes)
        self.numNodes = self.numBndNodes
        self.ndSize = self.numBndNodes
        self.nodes = np.array(boundaryNodes)
        
        
        self.numBndEdges = len(boundaryEdges)
        self.edgeNodes = np.array(boundaryEdges)
        self.numEdges = self.numBndEdges
        self.edSize = self.numBndEdges
        self.edgeElements = np.array([])
        self.edgeUnitNorms = np.array([])
        
        self.numTriEls = 0
        self.triElSize = 0
        self.triElements = np.array([])
        
        self.numQuadEls = 0
        self.quadElSize = 0
        self.quadElements = np.array([])
        
    ## !! check changes to createSweptMesh calls
    def createSweptMesh(self, sweepMethod, sweepElements, sweepDistance=1.0, point=[], axis=[], followNormal=False, destNodes=[], interpMethod='linear'):
        ## sweepMethod = inDirection, toPoint, fromPoint, toDestNodes, revolve
        """Object data modified: self.quadElements, self.nodes, self.quadElements
        Parameters
        ----------

        Returns
        -------
        nodes
        elements
        """ 
        nbNds = self.numBndNodes
        nbEds = self.numBndEdges
        try:
            totSweepEls = sum(sweepElements)
            ndSize = nbNds*(totSweepEls+1)
            stages = len(sweepElements)
            multiStage = True
        except:
            totSweepEls = sweepElements
            ndSize = nbNds*(sweepElements+1)
            stages = 1
            multiStage = False
        dimSpace = len(self.nodes[0])
        tmp = self.nodes.copy()
        self.nodes = np.zeros((ndSize,dimSpace))
        self.nodes[0:nbNds] = tmp
        self.ndSize = ndSize
        self.numNodes = nbNds
        
        if(self.numBndEdges == 0):
            n1 = np.array(range(nbNds-1))
            n2 = np.array(range(1,nbNds))
            self.edgeNodes = np.transpose(np.array([n1,n2]))
            self.numEdges = nbNds - 1
            self.numBndEdges = nbNds - 1
            self.edSize = nbNds - 1

        self.triElements = np.array([])
        self.triElSize = 0
        self.numTriEls = 0
        
        quadElSize = nbEds*totSweepEls
        self.quadElements = -np.ones((quadElSize,4),dtype=int)
        self.quadElSize = quadElSize
        self.numQuadEls = 0
        
        methString = 'inDirection toPoint fromPoint'
        if(sweepMethod in methString):
            ndDir = list()
            if(sweepMethod == 'inDirection'):
                mag = np.linalg.norm(axis)
                unitAxis = (1.0/mag)*np.array(axis)
                for i in range(0,self.numNodes):
                    ndDir.append(unitAxis)
            else:
                pAr = np.array(point)
                for i in range(0,self.numNodes):
                    if(sweepMethod == 'toPoint'):
                        vec = pAr - nd
                    else:
                        vec = nd - pAr
                    mag = np.linalg.norm(vec)
                    unitVec = (1.0/mag)*vec
                    ndDir.append(unitVec)
            rowNds = self.numNodes
            rowEls = self.numEdges
            stepLen = sweepDistance/sweepElements
            k = self.numNodes
            m = self.numQuadEls
            for i in range(0,sweepElements):
                for j in range(0,rowNds):
                    newNd = self.nodes[j] + (i+1)*stepLen*ndDir[j]
                    self.nodes[k] = newNd
                    k = k + 1
                for j in range(0,rowEls):
                    n1 = self.edgeNodes[j,0] + i*rowNds
                    n2 = self.edgeNodes[j,1] + i*rowNds
                    n3 = n2 + rowNds
                    n4 = n1 + rowNds
                    self.quadElements[m,:] = np.array([n1,n2,n3,n4])
                    m = m + 1
            self.numNodes = k
            self.numQuadEls = m
        elif(sweepMethod == 'toDestNodes'):
            nNds = self.numNodes
            nbNds = self.numBndNodes
            nEds = self.numEdges
            nQuad = self.numQuadEls
            if(not multiStage):
                sweepElements = [sweepElements]
                destNodes = [destNodes]
            if(interpMethod == 'linear'):
                prevDest = self.nodes.copy()
                for stg in range(0,stages):
                    dNds = np.array(destNodes[stg])
                    ndDir = list()
                    for ndi in range(0,nbNds):
                        vec = (1.0/sweepElements[stg])*(dNds[ndi] - prevDest[ndi])
                        ndDir.append(vec)
                    ndDir = np.array([ndDir])
                    for i in range(0,sweepElements[stg]):
                        for ndi in range(0,nbNds):
                            self.nodes[nNds] = self.nodes[ndi] + (i+1)*ndDir[ndi]
                            nNds = nNds + 1
                        for edi in range(0,self.numBndEdges):
                            n1 = self.edgeNodes[edi,0] + i*nbNds
                            n2 = self.edgeNodes[edi,1] + i*nbNds
                            n3 = n2 + nbNds
                            n4 = n1 + nbNds
                            self.quadElements[nQuad] = np.array([n1,n2,n3,n4])
                            nQuad = nQuad + 1
                    prevDest = dNds
            else:  ## Smooth interpolation
                xMat = np.zeros((nbNds,totSweepEls+1))
                yMat = np.zeros((nbNds,totSweepEls+1))
                pDest = (1.0/stages)*np.array(range(0,stages+1))
                pAll = (1.0/totSweepEls)*np.array(range(0,totSweepEls+1))
                for ndi in range(0,nbNds):
                    xDest = [self.nodes[ndi,0]]
                    yDest = [self.nodes[ndi,1]]
                    for dNds in destNodes:
                        xDest.append(dNds[ndi][0])
                        yDest.append(dNds[ndi][1])
                    xDest = np.array(xDest)
                    iFun = interpolate.interp1d(pDest,xDest,'cubic', axis=0,bounds_error=False,fill_value='extrapolate')
                    xAll = iFun(pAll)
                    xMat[ndi,:] = xAll
                    yDest = np.array(yDest)
                    iFun = interpolate.interp1d(pDest,yDest,'cubic', axis=0,bounds_error=False,fill_value='extrapolate')
                    yAll = iFun(pAll)
                    yMat[ndi,:] = yAll
                if(dimSpace == 3):
                    zMat = np.zeros((nbNds,totSweepEls+1))
                    for ndi in range(0,nbNds):
                        zDest = [self.nodes[ndi,2]]
                        for dNds in destNodes:
                            zDest.append(dNds[ndi][2])
                        zDest = np.array(zDest)
                        iFun = interpolate.interp1d(pDest,zDest,'cubic', axis=0,bounds_error=False,fill_value='extrapolate')
                        zAll = iFun(pAll)
                        zMat[ndi,:] = zAll
                for i in range(0,totSweepElements):
                    for ndi in range(0,nbNds):
                        if(dimSpace == 2):
                            self.nodes[nNds] = np.array([xMat[ndi,i+1],yMat[ndi,i+1]])
                        else:
                            self.nodes[nNds] = np.array([xMat[ndi,i+1],yMat[ndi,i+1],zMat[ndi,i+1]])
                        nNds = nNds + 1
                    for edi in range(0,self.numBndEdges):
                        n1 = self.edgeNodes[edi,0] + i*nbNds
                        n2 = self.edgeNodes[edi,1] + i*nbNds
                        n3 = n2 + nbNds
                        n4 = n1 + nbNds
                        self.quadElements[nQuad] = np.array([n1,n2,n3,n4])
                        nQuad = nQuad + 1
            self.numNodes = nNds
            self.numQuadEls = nQuad
        
        
        meshOut = dict()
        meshOut['nodes'] = self.nodes
        meshOut['elements'] = self.quadElements
        return meshOut
        
    def skewNodes(self):
        od = np.tan(np.pi/12.0)
        skewMat = np.array([[1.0,od],[od,1.0]])
        self.nodes = np.matmul(self.nodes,skewMat)
        
    def unskewNodes(self):
        od = np.tan(np.pi/12.0)
        skewMat = np.array([[1.0,od],[od,1.0]])
        invSkew = np.linalg.inv(skewMat)
        self.nodes = np.matmul(self.nodes,invSkew)
        
    def getBoundaryEdgeNormals(self):
        stepLen = self.minEdgeLen/np.sqrt(3.0)
        numSteps = int(np.ceil(self.edgeGL.xGSz*self.edgeGL.xRows/stepLen))
        yMin = self.edgeGL.yMin
        xMin = self.edgeGL.xMin
        xMarg = 0.6*self.maxEdgeLen
        for i in range(0,numSteps):
            xCrd = xMin + i*stepLen
            p1 = np.array([xCrd,yMin])
            v1 = np.array([0.0,1.0])
            Xns = list()
            nearEdges = self.edgeGL.findInXYMargin(p1,xMarg,-1)
            for edi in nearEdges:
                n1 = self.edgeNodes[edi,0]
                n2 = self.edgeNodes[edi,1]
                p2 = self.nodes[n1]
                v2 = self.nodes[n2] - p2
                Amat = np.array([[v1[0],-v2[0]],[v1[1],-v2[1]]])
                detA = np.linalg.det(Amat)
                if(detA != 0.0):
                    bVec = p2 - p1
                    soln = np.linalg.solve(Amat,bVec)
                    if(soln[1] > 0.0 and soln[1] < 1.0):
                        Xns.append([edi,soln[0]])
            iLen = len(Xns)
            for i in range(0,iLen-1):
                for j in range(0,iLen-1):
                    x1 = Xns[j]
                    x2 = Xns[j+1]
                    if(x2[1] < x1[1]):
                        Xns[j] = x2
                        Xns[j+1] = x1
            for i in range(0,iLen,2):
                edi = Xns[i][0]
                uN = self.edgeUnitNorms[edi]
                if(uN[1] < 0.0):
                    self.edgeUnitNorms[edi] = -uN
            for i in range(1,iLen,2):
                edi = Xns[i][0]
                uN = self.edgeUnitNorms[edi]
                if(uN[1] > 0.0):
                    self.edgeUnitNorms[edi] = -uN
                    
        numSteps = int(np.ceil(self.edgeGL.yGSz*self.edgeGL.yRows/stepLen))
        yMarg = 0.6*self.maxEdgeLen
        for i in range(0,numSteps):
            yCrd = yMin + i*stepLen
            p1 = np.array([xMin,yCrd])
            v1 = np.array([1.0,0.0])
            Xns = list()
            nearEdges = self.edgeGL.findInXYMargin(p1,-1,yMarg)
            for edi in nearEdges:
                n1 = self.edgeNodes[edi,0]
                n2 = self.edgeNodes[edi,1]
                p2 = self.nodes[n1]
                v2 = self.nodes[n2] - p2
                Amat = np.array([[v1[0],-v2[0]],[v1[1],-v2[1]]])
                detA = np.linalg.det(Amat)
                if(detA != 0.0):
                    bVec = p2 - p1
                    soln = np.linalg.solve(Amat,bVec)
                    if(soln[1] > 0.0 and soln[1] < 1.0):
                        Xns.append([edi,soln[0]])
            iLen = len(Xns)
            for i in range(0,iLen-1):
                for j in range(0,iLen-1):
                    x1 = Xns[j]
                    x2 = Xns[j+1]
                    if(x2[1] < x1[1]):
                        Xns[j] = x2
                        Xns[j+1] = x1
            for i in range(0,iLen,2):
                edi = Xns[i][0]
                uN = self.edgeUnitNorms[edi]
                if(uN[0] < 0.0):
                    self.edgeUnitNorms[edi] = -uN
            for i in range(1,iLen,2):
                edi = Xns[i][0]
                uN = self.edgeUnitNorms[edi]
                if(uN[0] > 0.0):
                    self.edgeUnitNorms[edi] = -uN        
    
    def unstructuredPrep(self,elType):
        nbNds = self.numBndNodes
        n_pi = nbNds/np.pi
        ndSize = int(4*n_pi*n_pi)
        tmp = self.nodes.copy()
        self.nodes = np.zeros((ndSize,2))
        self.nodes[0:nbNds,0:2] = tmp
        self.ndSize = ndSize
        self.numNodes = nbNds
        
        edSize = int(3*nbNds*n_pi)
        nbEds = self.numBndEdges
        if(nbEds == 0):
            n1 = np.array(range(nbNds))
            n2 = np.array(range(1,nbNds+1))
            n2[nbNds-1] = 0
            self.edgeNodes = -np.ones((edSize,2),dtype=int)
            self.edgeNodes[0:nbNds,0] = n1
            self.edgeNodes[0:nbNds,1] = n2
            nbEds = nbNds
            self.numBndEdges = nbEds
        else:
            tmp = self.edgeNodes.copy()
            self.edgeNodes = -np.ones((edSize,2),dtype=int)
            self.edgeNodes[0:nbEds] = tmp
        self.edgeElements = -np.ones((edSize,2),dtype=int)
        self.edgeUnitNorms = np.zeros((edSize,2))
        self.edSize = edSize
        self.numEdges = nbEds
                
        triElSize = int(2*nbNds*n_pi)
        self.triElements = -np.ones((triElSize,3),dtype=int)
        self.triElSize = triElSize
        self.numTriEls = 0
        
        quadElSize = int(nbNds*n_pi)
        self.quadElements = -np.ones((quadElSize,4),dtype=int)
        self.quadElSize = quadElSize
        self.numQuadEls = 0
        
        if(elType == 'quad'):
            self.skewNodes()
        minLen = 1.0e+100
        maxLen = 0.0
        avgLen = 0.0
        for edi in range(0,self.numEdges):
            n1 = self.edgeNodes[edi,0]
            n2 = self.edgeNodes[edi,1]
            vec = self.nodes[n1] - self.nodes[n2]
            ln = np.linalg.norm(vec)
            avgLen = avgLen + ln
            if(ln < minLen):
                minLen = ln
            if(ln > maxLen):
                maxLen = ln
            unitVec = (1.0/ln)*vec
            self.edgeUnitNorms[edi] = np.array([-unitVec[1],unitVec[0]])
        avgLen = avgLen/self.numEdges
        self.minEdgeLen = minLen
        self.maxEdgeLen = maxLen
        self.avgProjLen = 0.5*np.sqrt(3)*avgLen
        xMin = np.amin(self.nodes[0:self.numNodes,0])
        xMax = np.amax(self.nodes[0:self.numNodes,0])
        yMin = np.amin(self.nodes[0:self.numNodes,1])
        yMax = np.amax(self.nodes[0:self.numNodes,1])
        xLen = xMax - xMin
        yLen = yMax - yMin
        meshLen = np.sqrt(xLen*xLen + yLen*yLen)
        marg = 0.01*meshLen
        xMax = xMax + marg
        xMin = xMin - marg
        yMax = yMax + marg
        yMin = yMin - marg
        xLen = xMax - xMin
        yLen = yMax - yMin
        meshLen = np.sqrt(xLen*xLen + yLen*yLen)
        aL2 = avgLen*avgLen
        xS3 = xLen*aL2
        xGS = np.power(xS3,0.333333)
        yS3 = yLen*aL2
        yGS = np.power(yS3,0.333333)
        
        self.nodeGL = SpatialGridList2D(xMin,xMax,yMin,yMax,xGS,yGS)
        for ndi in range(0,self.numNodes):
            self.nodeGL.addEntry(ndi,self.nodes[ndi,:])
        
        self.edgeGL = SpatialGridList2D(xMin,xMax,yMin,yMax,xGS,yGS)
        for edi in range(0,self.numEdges):
            n1 = self.edgeNodes[edi,0]
            n2 = self.edgeNodes[edi,1]
            midPt = 0.5*(self.nodes[n1] + self.nodes[n2])
            self.edgeGL.addEntry(edi,midPt)
        
        self.triElGL = SpatialGridList2D(xMin,xMax,yMin,yMax,xGS,yGS)
        
        self.getBoundaryEdgeNormals()
        
    def edgesIntersect(self,e1Nds,e2Nds):
        e1n1 = e1Nds[0]
        e1n2 = e1Nds[1]
        e2n1 = e2Nds[0]
        e2n2 = e2Nds[1]
        p1 = self.nodes[e1n1]
        v1 = self.nodes[e1n2] - p1
        p2 = self.nodes[e2n1]
        v2 = self.nodes[e2n2] - p2
        Amat = np.array([[v1[0],-v2[0]],[v1[1],-v2[1]]])
        detA = np.linalg.det(Amat)
        if(detA != 0.0):
            bVec = p2 - p1
            soln = np.linalg.solve(Amat,bVec)
            if(soln[0] > 1.0e-6 and soln[0] < 0.999999 and soln[1] > 1.0e-6 and soln[1] < 0.999999):
                return True
            else:
                return False
        else:
            return False
            
    def ptInEl(self,pt,el):
        p1 = self.nodes[el[0]]
        v1 = self.nodes[el[1]] - p1
        v2 = self.nodes[el[2]] - p1
        Amat = np.array([[v1[0],v2[0]],[v1[1],v2[1]]])
        detA = np.linalg.det(Amat)
        if(detA != 0.0):
            bVec = pt - p1
            soln = np.linalg.solve(Amat,bVec)
            if(soln[0] > 1e-6 and soln[1] > 1e-6):
                solSum = soln[0] + soln[1]
                if(solSum < 0.999999):
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
            
    def violations(self,newEl):
        n1 = newEl[0]
        n2 = newEl[1]
        n3 = newEl[2]
        ed1 = np.array([n1,n2])
        ed2 = np.array([n2,n3])
        ed3 = np.array([n3,n1])
        cent = 0.333333*(self.nodes[n1] + self.nodes[n2] + self.nodes[n3])
        srchRad = 1.4*self.maxEdgeLen
        
        nearNds = self.nodeGL.findInRadius(cent,srchRad)
        for ndi in nearNds:
            inEl = self.ptInEl(self.nodes[ndi],newEl)
            if(inEl):
                return True
        
        nearEdges = self.edgeGL.findInRadius(cent,srchRad)
        for edi in nearEdges:
            intSect = self.edgesIntersect(self.edgeNodes[edi],ed1)
            if(intSect):
                return True
            intSect = self.edgesIntersect(self.edgeNodes[edi],ed2)
            if(intSect):
                return True
            intSect = self.edgesIntersect(self.edgeNodes[edi],ed3)
            if(intSect):
                return True
        
        srtNew = np.sort(newEl)
        nearEls = self.triElGL.findInRadius(cent,srchRad)
        for eli in nearEls:
            srtEi = np.sort(self.triElements[eli])
            if(all(srtNew == srtEi)):
                return True
            inEl = self.ptInEl(self.nodes[n1],self.triElements[eli])
            if(inEl):
                return True
            inEl = self.ptInEl(self.nodes[n2],self.triElements[eli])
            if(inEl):
                return True
            inEl = self.ptInEl(self.nodes[n3],self.triElements[eli])
            if(inEl):
                return True
                
        return False
        
    def createEdge(self,nds,el):
        midpt = 0.5*(self.nodes[nds[0]] + self.nodes[nds[1]])
        nearEdges = self.edgeGL.findInRadius(midpt,self.minEdgeLen)
        for nEi in nearEdges:    
            iNds = self.edgeNodes[nEi]
            if(nds[0] == iNds[0] and nds[1] == iNds[1]):
                self.edgeElements[nEi,1] = el
                return
            if(nds[0] == iNds[1] and nds[1] == iNds[0]):
                self.edgeElements[nEi,1] = el
                return
        k = self.numEdges
        self.edgeNodes[k] = nds
        self.edgeElements[k,0] = el
        edVec = self.nodes[nds[1]] - self.nodes[nds[0]]
        mag = np.linalg.norm(edVec)
        ## ---------------
        if(mag < 1.0e-12):
            outStr = 'nodes ' + str(nds)
            print(outStr)
            outStr = 'coords ' + str(self.nodes[nds[0]]) + '  ' + str(self.nodes[nds[1]])
            print(outStr)
        unitEV = (1.0/mag)*edVec
        unitNorm = np.array([-unitEV[1],unitEV[0]])
        for eNd in self.triElements[el]:
            if(eNd != nds[0] and eNd != nds[1]):
                vec = midpt - self.nodes[eNd]
                dp = np.dot(vec,unitNorm)
                if(dp > 0.0):
                    self.edgeUnitNorms[k] = unitNorm
                else:
                    self.edgeUnitNorms[k] = -unitNorm
        self.edgeGL.addEntry(k,midpt)
        self.numEdges = k + 1

    def adoptConnectedNode(self,edgeIndex,point,srchRad):
        eNds = self.edgeNodes[edgeIndex]
        midPt = 0.5*(self.nodes[eNds[0]] + self.nodes[eNds[1]])
        nearEdges = self.edgeGL.findInRadius(point,srchRad)
        for nEi in nearEdges:
            if(nEi != edgeIndex):
                eiNds = self.edgeNodes[nEi]
                commonNd = -1
                for i in range(0,2):
                    for j in range(0,2):
                        if(eNds[i] == eiNds[j]):
                            commonNd = j
                if(commonNd != -1):
                    if(commonNd == 0):
                        commonNdi = eiNds[0]
                        unComNdi = eiNds[1]
                    else:
                        commonNdi = eiNds[1]
                        unComNdi = eiNds[0]
                    unComNdPt = self.nodes[unComNdi]
                    vec = point - unComNdPt
                    dist = np.linalg.norm(vec)
                    if(dist < srchRad):
                        vec = unComNdPt - midPt
                        dp = np.dot(vec,self.edgeUnitNorms[edgeIndex])
                        if(dp > 0.0):
                            newEl = np.array([eNds[0],eNds[1],unComNdi])
                            viol = self.violations(newEl)
                            if(not viol):
                                k = self.numTriEls
                                self.triElements[k] = newEl
                                cent = 0.33333333*(self.nodes[newEl[0]] + self.nodes[newEl[1]] + self.nodes[newEl[2]])
                                self.triElGL.addEntry(k,cent)
                                self.edgeElements[edgeIndex,1] = k
                                self.edgeElements[nEi,1] = k
                                if(eNds[0] == commonNdi):
                                    newNds = np.array([eNds[1],unComNdi])
                                else:
                                    newNds = np.array([eNds[0],unComNdi])
                                self.createEdge(newNds,k)
                                self.numTriEls = k + 1
                                return True
        return False
                
    def adoptAnyNode(self,edgeIndex,point,srchRad):
        eNds = self.edgeNodes[edgeIndex]
        midPt = 0.5*(self.nodes[eNds[0]] + self.nodes[eNds[1]])
        nearNds = self.nodeGL.findInRadius(point,srchRad)
        for ndi in nearNds:
            if(ndi not in eNds):
                vec = self.nodes[ndi] - point
                dist = np.linalg.norm(vec)
                if(dist < srchRad):
                    vec = self.nodes[ndi] - midPt
                    dp = np.dot(vec,self.edgeUnitNorms[edgeIndex])
                    if(dp > 0.0):
                        newEl = np.array([eNds[0],eNds[1],ndi])
                        viol = self.violations(newEl)
                        if(not viol):
                            k = self.numTriEls
                            self.triElements[k] = newEl
                            cent = 0.33333333*(self.nodes[newEl[0]] + self.nodes[newEl[1]] + self.nodes[newEl[2]])
                            self.triElGL.addEntry(k,cent)
                            self.edgeElements[edgeIndex,1] = k
                            newNds = np.array([eNds[0],ndi])
                            self.createEdge(newNds,k)
                            newNds = np.array([eNds[1],ndi])
                            self.createEdge(newNds,k)
                            self.numTriEls = k + 1
                            return True
        return False
        
    def createNode(self,edgeIndex,point):
        eNds = self.edgeNodes[edgeIndex]
        n = self.numNodes
        self.nodes[n] = point
        newEl = np.array([eNds[0],eNds[1],n])
        viol = self.violations(newEl)
        if(not viol):
            k = self.numTriEls
            self.triElements[k] = newEl
            cent = 0.33333333*(self.nodes[newEl[0]] + self.nodes[newEl[1]] + self.nodes[newEl[2]])
            self.triElGL.addEntry(k,cent)
            self.edgeElements[edgeIndex,1] = k
            newNds = np.array([eNds[0],n])
            self.createEdge(newNds,k)
            newNds = np.array([eNds[1],n])
            self.createEdge(newNds,k)
            self.nodeGL.addEntry(n,point)
            self.numTriEls = k + 1
            self.numNodes = n + 1
            return True
        else:
            return False
            
    def distributeNodes(self):
        dim = 2*self.numNodes
        Dmat = np.zeros(dim)
        bDim = 2*self.numBndNodes
        Dmat[0:bDim] = 100000.0
        Pmat = 10.0*np.ones(dim) + Dmat
        Pinv = np.zeros(dim)
        Pinv[0:bDim] = 9.999e-6
        Pinv[bDim:dim] = 0.1
        rhs = np.zeros(dim)
        j = 0
        for bni in range(0,self.numBndNodes):
            rhs[j] = 100000.0*self.nodes[bni,0]
            rhs[j+1] = 100000.0*self.nodes[bni,1]
            j = j + 2
        
        nEls = self.numTriEls
        elWt = np.zeros(nEls)
        for eli in range(0,nEls):
            ni = self.triElements[eli]
            v1 = self.nodes[ni[1]] - self.nodes[ni[0]]
            v2 = self.nodes[ni[2]] - self.nodes[ni[0]]
            cp = v1[0] * v2[1] - v1[1] * v2[0]
            elWt[eli] = np.abs(cp)        
        avgWt = np.mean(elWt)
        elWt = (1.0/avgWt)*elWt
        
        elMat = np.zeros((6,6))
        elMat[0,0] = 2.0
        elMat[0,2] = -1.0
        elMat[0,4] = -1.0
        for i in range(1,6):
            elMat[i,i:6] = elMat[0,0:6-i]
        for i in range(0,5):
            elMat[i+1:6,i] = elMat[i,i+1:6]
        
        xVec = np.zeros(dim)
        gVec = -rhs
        wVec = np.multiply(Pinv,gVec)
        hVec = -wVec
        zVec = np.zeros(dim)
        res = np.dot(gVec,wVec)
        i = 0
        while(res > 1e-12 and i < dim):
            zVec[:] = 0.0
            for eli in range(0,nEls):
                inT2 = 2*self.triElements[eli]
                vecInd = [inT2[0],inT2[0]+1,inT2[1],inT2[1]+1,inT2[2],inT2[2]+1]
                elH = hVec[vecInd]
                elMati = elWt[eli]*elMat
                elZ = np.matmul(elMati,elH)
                zVec[vecInd] = zVec[vecInd] + elZ
            zVec = zVec + np.multiply(Dmat,hVec)
            alpha = res/np.dot(hVec,zVec)
            xVec = xVec + alpha*hVec
            gVec = gVec + alpha*zVec
            wVec = np.multiply(Pinv,gVec)
            rNext = np.dot(gVec,wVec)
            beta = rNext/res
            res = rNext
            hVec = -wVec
            i = i + 1
        
        for i in range(0,self.numNodes):
            j = i*2
            self.nodes[i] = xVec[j:j+2]
            
    def mergePairsAbove(self,edgeFactor,elElim,elLongEdge):
        nQuad = self.numQuadEls
        for edi in range(0,self.numEdges):
            n1 = self.edgeNodes[edi,0]
            n2 = self.edgeNodes[edi,1]
            edLen = np.linalg.norm(self.nodes[n1]-self.nodes[n2])
            el1 = self.edgeElements[edi,0]
            el2 = self.edgeElements[edi,1]
            if(el1 != -1 and el2 != -1):
                if(elElim[el1] == 0 and elElim[el2] == 0):
                    if(edLen > edgeFactor*elLongEdge[el1] and edLen > edgeFactor*elLongEdge[el2]):
                        elElim[el1] = 1
                        elElim[el2] = 1
                        quadNds = np.array([n1,-1,n2,-1])
                        e1Nds = self.triElements[el1]
                        e2Nds = self.triElements[el2]
                        for i in range(0,3):
                            if(e1Nds[i] != n1 and e1Nds[i] != n2):
                                quadNds[1] = e1Nds[i]
                            if(e2Nds[i] != n1 and e2Nds[i] != n2):
                                quadNds[3] = e2Nds[i]
                        self.quadElements[nQuad] = quadNds
                        nQuad = nQuad + 1
        self.numQuadEls = nQuad
        return elElim
            
    def mergeTriEls(self,elType):
        nNds = self.numNodes
        nbNds = self.numBndNodes
        nEls = self.numTriEls
        nQuad = self.numQuadEls
        elElim = np.zeros(self.triElSize,dtype=int)
        ndElems = -np.ones((nNds,6),dtype=int)
        ndElems[:,0] = 0
        ndElim = np.zeros(nNds,dtype=int)
        
        for eli in range(0,nEls):
            for nd in self.triElements[eli]:
                j = ndElems[nd,0]
                if(j < 5):
                    j = j + 1
                    ndElems[nd,j] = eli
                    ndElems[nd,0] = j
        
        for ndi in range(nbNds,nNds):
            if(ndElems[ndi,4] == -1):
                abrt = False
                for el in ndElems[ndi,1:4]:
                    if(elElim[el] == 1):
                        abrt = True
                if(not abrt):    
                    ndElim[ndi] = 1
                    newElNds = list()
                    for el in ndElems[ndi,1:4]:
                        elElim[el] = 1
                        newElNds.extend(self.triElements[el])
                    srtedNds = np.sort(newElNds)
                    finalNds = list()
                    for i in range(0,8):
                        j = srtedNds[i]
                        if(j != ndi and srtedNds[i+1] == j):
                            finalNds.append(j)
                    self.triElements[nEls] = np.array(finalNds,dtype=int)
                    nEls = nEls + 1
            elif(ndElems[ndi,5] == -1):
                abrt = False
                for el in ndElems[ndi,1:5]:
                    if(elElim[el] == 1):
                        abrt = True
                if(not abrt):   
                    ndElim[ndi] = 1
                    newElNds = list()
                    for el in ndElems[ndi,1:4]:
                        elElim[el] = 1
                        newElNds.extend(self.triElements[el])
                    elElim[ndElems[ndi,4]] = 1
                    srtedNds = np.sort(newElNds)
                    nds12 = list()
                    for i in range(0,8):
                        j = srtedNds[i]
                        if(j != ndi and j == srtedNds[i+1]):
                            nds12.append(j)
                    nds34 = list()
                    for i in range(0,8):
                        j = srtedNds[i]
                        if(j != ndi and j not in nds12):
                            nds34.append(j)
                    n1 = nds12[0]
                    n2 = nds12[1]
                    n3 = nds34[0]
                    n4 = nds34[1]
                    v1 = self.nodes[n2] - self.nodes[n1]
                    v2 = self.nodes[n3] - self.nodes[n4]
                    dp = np.dot(v1,v2)
                    if(dp > 0.0):
                        if(elType == 'quad'):
                            self.quadElements[nQuad] = np.array([n1,n2,n3,n4])
                            nQuad = nQuad + 1
                        else:
                            self.triElements[nEls] = np.array([n1,n2,n3])
                            nEls = nEls + 1
                            self.triElements[nEls] = np.array([n1,n3,n4])
                            nEls = nEls + 1
                    else:
                        if(elType == 'quad'):
                            self.quadElements[nQuad] = np.array([n1,n2,n4,n3])
                            nQuad = nQuad + 1
                        else:
                            self.triElements[nEls] = np.array([n1,n2,n4])
                            nEls = nEls + 1
                            self.triElements[nEls] = np.array([n1,n4,n3])
                            nEls = nEls + 1
        self.numTriEls = nEls
        self.numQuadEls = nQuad
        
        if(elType == 'quad'):
            elLongEdge = np.zeros(nEls)
        
            for eli in range(0,nEls):
                nds = self.triElements[eli]
                longEd = np.linalg.norm(self.nodes[nds[0]] - self.nodes[nds[1]])
                eLen = np.linalg.norm(self.nodes[nds[1]] - self.nodes[nds[2]])
                if(eLen > longEd):
                    longEd = eLen
                eLen = np.linalg.norm(self.nodes[nds[2]] - self.nodes[nds[0]])
                if(eLen > longEd):
                    longEd = eLen
                elLongEdge[eli] = longEd
        
            ## Merge quad-forming pairs of triangles
            
            elElim = self.mergePairsAbove(0.99,elElim,elLongEdge)
            elElim = self.mergePairsAbove(0.85,elElim,elLongEdge)
            elElim = self.mergePairsAbove(0.75,elElim,elLongEdge)
            nQuad = self.numQuadEls
        
        finalNodes = list()
        ndNewInd = -np.ones(nNds,dtype=int)
        ndi = 0
        for ni in range(0,nNds):
            if(ndElim[ni] == 0):
                finalNodes.append(self.nodes[ni])
                ndNewInd[ni] = ndi
                ndi = ndi + 1
        self.nodes = np.array(finalNodes)
        self.numNodes = len(self.nodes)
        self.ndSize = self.numNodes
        
        for eli in range(0,nEls):
            if(elElim[eli] == 0):
                for j in range(0,3):
                    nd = self.triElements[eli,j]
                    self.triElements[eli,j] = ndNewInd[nd]
        
        for eli in range(0,nQuad):
            for j in range(0,4):
                nd = self.quadElements[eli,j]
                self.quadElements[eli,j] = ndNewInd[nd]
        
        newTEind = list()
        for i in range(0,nEls):
            if(elElim[i] == 0):
                newTEind.append(i)

        self.triElements = self.triElements[newTEind]
        self.numTriEls = len(newTEind)
        self.triElSize = self.numTriEls
            
    def unstructuredPost(self,elType):
        #self.plot2DMesh()
        self.distributeNodes()
        #self.plot2DMesh()
        if(elType == 'quad'):
            self.unskewNodes()
        #self.plot2DMesh()
        self.mergeTriEls(elType)
        #self.plot2DMesh()
                        
    ## !! rename any calls to creating unstructured mesh as necessary, createPlanarMesh
    def createUnstructuredMesh(self,elType):
        self.unstructuredPrep(elType)
        
        elsCreated = True
        while(elsCreated):
            # if(self.numTriEls > 0):
                # self.plot2DMesh()
            # cnt = input('continue?\n')
            # if(cnt != 'y'):
                # break
            elsCreated = False
            nEd = self.numEdges
            for edi in range(0,nEd):
                if(self.edgeElements[edi,1] == -1):
                    n1 = self.edgeNodes[edi,0]
                    n2 = self.edgeNodes[edi,1]
                    uNorm = self.edgeUnitNorms[edi]
                    midPt = 0.5*(self.nodes[n1] + self.nodes[n2])
                    vec = self.nodes[n1] - self.nodes[n2]
                    edLen = np.linalg.norm(vec)
                    projLen = 0.2165*edLen + 0.75*self.avgProjLen ## 0.25(sqrt(3)/2)edLen + 0.75avgProjLen
                    srchPt = midPt + 0.5*projLen*uNorm
                    srchRad = 0.5*np.sqrt(edLen*edLen + projLen*projLen)
                    found = self.adoptConnectedNode(edi,srchPt,srchRad)
                    if(not found):
                        srchPt = midPt + projLen*uNorm
                        found = self.adoptConnectedNode(edi,srchPt,srchRad)
                    if(not found):
                        srchPt = midPt + 0.5*projLen*uNorm
                        found = self.adoptAnyNode(edi,srchPt,srchRad)
                    if(not found):
                        srchPt = midPt + projLen*uNorm
                        found = self.adoptAnyNode(edi,srchPt,srchRad)
                    if(not found):
                        srchPt = midPt + projLen*uNorm
                        self.createNode(edi,srchPt)
                    if(not found):
                        srchPt = midPt + 0.5*projLen*uNorm
                        self.createNode(edi,srchPt)
                    if(found):
                        elsCreated = True
            for edi in range(0,self.numEdges):
                uNMag = np.linalg.norm(self.edgeUnitNorms[edi])
                if(uNMag < 1.0e-6):
                    n1 = self.edgeNodes[edi,0]
                    n2 = self.edgeNodes[edi,1]
                    midPt = 0.5*(self.nodes[n1] + self.nodes[n2])
                    vec = self.nodes[n2] - self.nodes[n1]
                    mag = np.linalg.norm(vec)
                    unitVec = (1.0/mag)*vec
                    uNorm = np.array([-unitVec[1],unitVec[0]])
                    eL1 = self.edgeElements[edi,0]
                    for i in range(0,3):
                        ni = self.triElements[eL1,i]
                        if(ni != n1 and ni != n2):
                            vec = midPt - self.nodes[ni]
                            dp = np.dot(vec,uNorm)
                            if(dp > 0.0):
                                self.edgeUnitNorms[edi] = uNorm
                            else:
                                self.edgeUnitNorms[edi] = -uNorm
        
        self.unstructuredPost(elType)
        
        meshOut = dict()
        meshOut['nodes'] = self.nodes
        totalEls = self.numTriEls + self.numQuadEls
        allEls = -np.ones((totalEls,4),dtype=int)
        allEls[0:self.numTriEls,0:3] = self.triElements[0:self.numTriEls]
        allEls[self.numTriEls:totalEls,:] = self.quadElements[0:self.numQuadEls]
        meshOut['elements'] = allEls
        
        return meshOut
        
    def plot2DMesh(self):
        xLst = self.nodes[0:self.numNodes,0]
        yLst = self.nodes[0:self.numNodes,1]
        zLst = np.zeros(self.numNodes)
        value = list()
        v1 = list()
        v2 = list()
        v3 = list()
        for i in range(0,self.numTriEls):
            v1.append(self.triElements[i,0])
            v2.append(self.triElements[i,1])
            v3.append(self.triElements[i,2])
            value.append(np.sin(i))
        for i in range(0,self.numQuadEls):
            v1.append(self.quadElements[i,0])
            v2.append(self.quadElements[i,1])
            v3.append(self.quadElements[i,2])
            value.append(np.sin(i))
            v1.append(self.quadElements[i,0])
            v2.append(self.quadElements[i,2])
            v3.append(self.quadElements[i,3])
            value.append(np.sin(i))
        fig = go.Figure(data=[
            go.Mesh3d(
                x=xLst,
                y=yLst,
                z=zLst,
                colorbar_title = '',
                colorscale=[[0.0, 'blue'],
                            [0.5, 'yellow'],
                            [1.0, 'red']],
                intensity=value,
                intensitymode='cell',
                i=v1,
                j=v2,
                k=v3,
                name='',
                showscale=True
            )
        ])

        fig.show()