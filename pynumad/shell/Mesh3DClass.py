import numpy as np
from scipy import interpolate
from pynumad.shell.SpatialGridList3DClass import *

class Mesh3D():

    def __init__(self,boundaryNodes,boundaryFaces=[]):
        self.nodeGL = None
        self.faceGL = None
        self.tetElGL = None
        
        self.minFaceArea = 0.0
        self.maxFaceArea = 1.0
        self.avgProjLen = 1.0
        
        self.numBndNodes = len(boundaryNodes)
        self.numNodes = self.numBndNodes
        self.ndSize = self.numBndNodes
        self.nodes = np.array(boundaryNodes)
        
        self.numBndFaces = len(boundaryFaces)
        self.faceNodes = np.array(boundaryFaces)
        self.numFaces = self.numBndFaces
        self.faceSize = self.numBndFaces
        self.faceElements = np.array([])
        self.faceUnitNorms = np.array([])
        
        self.numTetEls = 0
        self.tetElSize = 0
        self.tetElements = np.array([])
        
        self.numWedgeEls = 0
        self.wedgeElSize = 0
        self.wedgeElements = np.array([])
        
        self.numHexEls = 0
        self.hexElSize = 0
        self.hexElements = np.array([])
        
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
        tmp = self.nodes.copy()
        self.nodes = np.zeros((ndSize,3))
        self.nodes[0:nbNds] = tmp
        self.ndSize = ndSize
        self.numNodes = nbNds
        
        nbFcs = self.numBndFaces
        elSize = nbFcs*totSweepEls
        
        self.wedgeElements = -np.ones((elSize,6),dtype=int)
        self.wedgeElSize = elSize
        self.numWedgeEls = 0
        
        self.hexElements = -np.ones((elSize,8),dtype=int)
        self.hexElSize = elSize
        self.numHexEls = 0
        
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
            rowEls = self.numFaces
            stepLen = sweepDistance/sweepElements
            nNds = self.numNodes
            wE = self.numWedgeEls
            hE = self.numHexEls
            for i in range(0,sweepElements):
                for j in range(0,rowNds):
                    newNd = self.nodes[j] + (i+1)*stepLen*ndDir[j]
                    self.nodes[nNds] = newNd
                    nNds = nNds + 1
                for j in range(0,rowEls):
                    n1 = self.faceNodes[j,0] + i*rowNds
                    n2 = self.faceNodes[j,1] + i*rowNds
                    n3 = self.faceNodes[j,2] + i*rowNds
                    if(self.faceNodes[j,3] == -1):
                        n4 = n1 + rowNds
                        n5 = n2 + rowNds
                        n6 = n3 + rowNds
                        self.wedgeElements[wE] = np.array([n1,n2,n3,n4,n5,n6])
                        wE = wE + 1
                    else:
                        n4 = self.faceNodes[j,3] + i*rowNds
                        n5 = n1 + rowNds
                        n6 = n2 + rowNds
                        n7 = n3 + rowNds
                        n8 = n4 + rowNds
                        self.hexElements[hE] = np.array([n1,n2,n3,n4,n5,n6,n7,n8])
                        hE = hE + 1
            
            self.numNodes = nNds
            self.numWedgeEls = wE
            self.numHexEls = hE
        elif(sweepMethod == 'toDestNodes'):
            nNds = self.numNodes
            nbNds = self.numBndNodes
            wE = self.numWedgeEls
            hE = self.numHexEls
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
                        for fci in range(0,self.numBndFaces):
                            n1 = self.faceNodes[fci,0] + i*nbNds
                            n2 = self.faceNodes[fci,1] + i*nbNds
                            n3 = self.faceNodes[fci,2] + i*nbNds
                            if(self.faceNodes[fci,3] == -1):
                                n4 = n1 + nbNds
                                n5 = n2 + nbNds
                                n6 = n3 + nbNds
                                self.wedgeElements[wE] = np.array([n1,n2,n3,n4,n5,n6])
                                wE = wE + 1
                            else:
                                n4 = self.faceNodes[fci,3] + i*nbNds
                                n5 = n1 + nbNds
                                n6 = n2 + nbNds
                                n7 = n3 + nbNds
                                n8 = n4 + nbNds
                                self.hexElements[hE] = np.array([n1,n2,n3,n4,n5,n6,n7,n8])
                                hE = hE + 1
                    prevDest = dNds
            else:  ## Smooth interpolation
                xMat = np.zeros((nbNds,totSweepEls+1))
                yMat = np.zeros((nbNds,totSweepEls+1))
                zMat = np.zeros((nbNds,totSweepEls+1))
                pDest = (1.0/stages)*np.array(range(0,stages+1))
                pAll = (1.0/totSweepEls)*np.array(range(0,totSweepEls+1))
                for ndi in range(0,nbNds):
                    xDest = [self.nodes[ndi,0]]
                    yDest = [self.nodes[ndi,1]]
                    zDest = [self.nodes[ndi,2]]
                    for dNds in destNodes:
                        xDest.append(dNds[ndi][0])
                        yDest.append(dNds[ndi][1])
                        zDest.append(dNds[ndi][2])
                    xDest = np.array(xDest)
                    iFun = interpolate.interp1d(pDest,xDest,'cubic', axis=0,bounds_error=False,fill_value='extrapolate')
                    xAll = iFun(pAll)
                    xMat[ndi,:] = xAll
                    yDest = np.array(yDest)
                    iFun = interpolate.interp1d(pDest,yDest,'cubic', axis=0,bounds_error=False,fill_value='extrapolate')
                    yAll = iFun(pAll)
                    yMat[ndi,:] = yAll
                    zDest = np.array(zDest)
                    iFun = interpolate.interp1d(pDest,zDest,'cubic', axis=0,bounds_error=False,fill_value='extrapolate')
                    zAll = iFun(pAll)
                    zMat[ndi,:] = zAll
                for i in range(0,totSweepEls):
                    for ndi in range(0,nbNds):
                        self.nodes[nNds] = np.array([xMat[ndi,i+1],yMat[ndi,i+1],zMat[ndi,i+1]])
                        nNds = nNds + 1
                    for fci in range(0,self.numBndFaces):
                        n1 = self.faceNodes[fci,0] + i*nbNds
                        n2 = self.faceNodes[fci,1] + i*nbNds
                        n3 = self.faceNodes[fci,2] + i*nbNds
                        if(self.faceNodes[fci,3] == -1):
                            n4 = n1 + nbNds
                            n5 = n2 + nbNds
                            n6 = n3 + nbNds
                            self.wedgeElements[wE] = np.array([n1,n2,n3,n4,n5,n6])
                            wE = wE + 1
                        else:
                            n4 = self.faceNodes[fci,3] + i*nbNds
                            n5 = n1 + nbNds
                            n6 = n2 + nbNds
                            n7 = n3 + nbNds
                            n8 = n4 + nbNds
                            self.hexElements[hE] = np.array([n1,n2,n3,n4,n5,n6,n7,n8])
                            hE = hE + 1
            self.numNodes = nNds
            self.numWedgeEls = wE
            self.numHexEls = hE
        
        for eli in range(0,self.numWedgeEls):
            n1 = self.wedgeElements[eli,0]
            n2 = self.wedgeElements[eli,1]
            n3 = self.wedgeElements[eli,2]
            n4 = self.wedgeElements[eli,3]
            v1 = self.nodes[n2] - self.nodes[n1]
            v2 = self.nodes[n3] - self.nodes[n1]
            v3 = self.nodes[n4] - self.nodes[n1]
            mat = np.array([v1,v2,v3])
            detM = np.linalg.det(mat)
            if(detM < 0.0):
                self.wedgeElements[eli,1] = n3
                self.wedgeElements[eli,2] = n2
                sw = self.wedgeElements[eli,4]
                self.wedgeElements[eli,4] = self.wedgeElements[eli,5]
                self.wedgeElements[eli,5] = sw
                
        for eli in range(0,self.numHexEls):
            n1 = self.hexElements[eli,0]
            n2 = self.hexElements[eli,1]
            n4 = self.hexElements[eli,3]
            n5 = self.hexElements[eli,4]
            v1 = self.nodes[n2] - self.nodes[n1]
            v2 = self.nodes[n4] - self.nodes[n1]
            v3 = self.nodes[n5] - self.nodes[n1]
            mat = np.array([v1,v2,v3])
            detM = np.linalg.det(mat)
            if(detM < 0.0):
                self.hexElements[eli,1] = n4
                self.hexElements[eli,3] = n2
                sw = self.hexElements[eli,5]
                self.hexElements[eli,5] = self.hexElements[eli,7]
                self.hexElements[eli,7] = sw
                
        
        meshOut = dict()
        meshOut['nodes'] = self.nodes
        totEls = self.numWedgeEls + self.numHexEls
        allEls = -np.ones((totEls,8),dtype=int)
        if(self.numWedgeEls > 0):
            allEls[0:self.numWedgeEls,0:6] = self.wedgeElements[0:self.numWedgeEls]
        if(self.numHexEls > 0):
            allEls[self.numWedgeEls:totEls,0:8] = self.hexElements[0:self.numHexEls]
        meshOut['elements'] = allEls
        return meshOut