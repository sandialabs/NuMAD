from scipy import interpolate
from pynumad.shell.Segment2DClass import *
from pynumad.shell.Boundary2DClass import *
from pynumad.shell.Mesh2DClass import *
import pynumad.shell.MeshTools as mt

class ShellRegion: 
    """
    Attributes
    -----------
    type : str
    keyPts : list
    edgeEls : list
    """
    def __init__(self, regType, keyPoints, numEdgeEls, natSpaceCrd=[], elType='quad', meshMethod='free'): 

        self.regType = regType
        self.keyPts = np.array(keyPoints)
        self.edgeEls = numEdgeEls
        if(len(natSpaceCrd) == 0):
            if(regType == 'quad1'):
                self.natSpaceCrd = np.array([[-1.0,-1.0],[1.0,-1.0],[1.0,1.0],[-1.0,1.0]])
            elif(regType == 'quad2'):
                self.natSpaceCrd = np.array([[-1.0,-1.0],[1.0,-1.0],[1.0,1.0],[-1.0,1.0],
                    [0.0,-1.0],[1.0,0.0],[0.0,1.0],[-1.0,0.0],[0.0,0.0]])
            elif(regType == 'quad3'):
                r3 = 1.0/3.0
                self.natSpaceCrd = np.array([[-1.0,-1.0],[1.0,-1.0],[1.0,1.0],[-1.0,1.0],
                    [-r3,-1.0],[r3,-1.0],[1.0,-r3],[1.0,r3],[r3,1.0],[-r3,1.0],[-1.0,r3],[-1.0,-r3],
                    [-r3,-r3],[r3,-r3],[r3,r3],[-r3,r3]])
        else:
            self.natSpaceCrd = np.array(natSpaceCrd)
        self.elType = elType
        self.meshMethod = meshMethod
    
    def createShellMesh(self):
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
        if(self.meshMethod == 'structured'):
            if('quad' in self.regType):
                ee = self.edgeEls
                if(ee[0] >= ee[2]):
                    xNodes = ee[0] + 1
                else:
                    xNodes = ee[2] + 1
                if(ee[1] >= ee[3]):
                    yNodes = ee[1] + 1
                else:
                    yNodes = ee[3] + 1
                totNds = xNodes*yNodes
                seg = Segment2D('line',[[-1.0,-1.0],[1.0,-1.0]],(xNodes-1))
                bnd = seg.getNodesEdges()
                mesh = Mesh2D(bnd['nodes'],bnd['edges'])
                mData = mesh.createSweptMesh('inDirection',(yNodes-1),sweepDistance=2.0,axis=[0.0,1.0])

                moved = False
                if(self.edgeEls[0] < self.edgeEls[2]):
                    seg = Segment2D('line',[[-1.0,-1.0],[1.0,-1.0]],self.edgeEls[0])
                    bnd = seg.getNodesEdges()
                    segNds = bnd['nodes']
                    meshNds = mData['nodes']
                    for ndi in range(0,xNodes):
                        minDist = 2.0
                        for sN in segNds:
                            vec = meshNds[ndi] - sN
                            dist = np.linalg.norm(vec)
                            if(dist < minDist):
                                minDist = dist
                                minPt = sN
                        meshNds[ndi] = minPt
                    mData['nodes'] = meshNds
                    moved = True
                elif(self.edgeEls[2] < self.edgeEls[0]):
                    seg = Segment2D('line',[[-1.0,1.0],[1.0,1.0]],self.edgeEls[2])
                    bnd = seg.getNodesEdges()
                    segNds = bnd['nodes']
                    meshNds = mData['nodes']
                    for ndi in range((totNds-xNodes),totNds):
                        minDist = 2.0
                        for sN in segNds:
                            vec = meshNds[ndi] - sN
                            dist = np.linalg.norm(vec)
                            if(dist < minDist):
                                minDist = dist
                                minPt = sN
                        meshNds[ndi] = minPt
                    mData['nodes'] = meshNds
                    moved = True
                if(self.edgeEls[1] < self.edgeEls[3]):
                    seg = Segment2D('line',[[1.0,-1.0],[1.0,1.0]],self.edgeEls[1])
                    bnd = seg.getNodesEdges()
                    segNds = bnd['nodes']
                    meshNds = mData['nodes']
                    for ndi in range((xNodes-1),totNds,xNodes):
                        minDist = 2.0
                        for sN in segNds:
                            vec = meshNds[ndi] - sN
                            dist = np.linalg.norm(vec)
                            if(dist < minDist):
                                minDist = dist
                                minPt = sN
                        meshNds[ndi] = minPt
                    mData['nodes'] = meshNds
                    moved = True     
                elif(self.edgeEls[3] < self.edgeEls[1]):
                    seg = Segment2D('line',[[-1.0,1.0],[-1.0,-1.0]],self.edgeEls[3])
                    bnd = seg.getNodesEdges()
                    segNds = bnd['nodes']
                    meshNds = mData['nodes']
                    for ndi in range(0,totNds,xNodes):
                        minDist = 2.0
                        for sN in segNds:
                            vec = meshNds[ndi] - sN
                            dist = np.linalg.norm(vec)
                            if(dist < minDist):
                                minDist = dist
                                minPt = sN
                        meshNds[ndi] = minPt
                    mData['nodes'] = meshNds
                    moved = True
                
                if(moved):                    
                    mData = mt.mergeDuplicateNodes(mData)
                    elLst = mData['elements']
                    ndLst = mData['nodes']     
                    for eli in range(0,len(elLst)):
                        srted = np.sort(elLst[eli])
                        for i in range(0,3):
                            if(srted[i+1] == srted[i]):
                                srted[i+1] = srted[3]
                                srted[3] = -1
                                elLst[eli] = srted
                        if(elLst[eli,3] == -1):
                            n1 = elLst[eli,0]
                            n2 = elLst[eli,1]
                            n3 = elLst[eli,2]
                            v1 = ndLst[n2] - ndLst[n1]
                            v2 = ndLst[n3] - ndLst[n1]
                            k = v1[0]*v2[1] - v1[1]*v2[0]
                            if(k < 0.0):
                                elLst[eli,1] = n3
                                elLst[eli,2] = n2
                else:
                    elLst = mData['elements']
                    ndLst = mData['nodes'] 

                XYZ = self.XYZCoord(ndLst)
                
                mData['nodes'] = XYZ
                mData['elements'] = elLst
                return mData

            else:
                raise Exception('Only quadrilateral shell regions can use the structured meshing option')

        else:
            bndData = self.initialBoundary()
            mesh = Mesh2D(bndData['nodes'],bndData['elements'])
            mData = mesh.createUnstructuredMesh(self.elType)
            XYZ = self.XYZCoord(mData['nodes'])
            mData['nodes'] = XYZ
            return mData
        
    def initialBoundary(self):
        """ Object data modified: none
        Parameters
        ----------

        Returns
        -------
        nodes
        edges
        """
        if 'quad' in self.regType:
            bnd = Boundary2D()
            bnd.addSegment('line',[[-1.0,-1.0],[1.0,-1.0]],self.edgeEls[0])
            bnd.addSegment('line',[[1.0,-1.0],[1.0,1.0]],self.edgeEls[1])
            bnd.addSegment('line',[[1.0,1.0],[-1.0,1.0]],self.edgeEls[2])
            bnd.addSegment('line',[[-1.0,1.0],[-1.0,-1.0]],self.edgeEls[3])
            bData = bnd.getBoundaryMesh()
            return bData
        elif 'tri' in self.regType:
            bnd = Boundary2D()
            bnd.addSegment('line',[[0.0,0.0],[1.0,0.0]],self.edgeEls[0])
            bnd.addSegment('line',[[1.0,0.0],[0.0,1.0]],self.edgeEls[1])
            bnd.addSegment('line',[[0.0,1.0],[0.0,0.0]],self.edgeEls[2])
            bData = bnd.getBoundaryMesh()
            return bData
        elif 'sphere' in self.regType:
            pi_2 = 0.5*np.pi
            bnd = Boundary2D()
            bnd.addSegment('arc',[[pi_2,0.0],[-pi_2,0.0],[pi_2,0.0]],self.edgeEls[0])
            bData = bnd.getBoundaryMesh()
            return bData
        
        
    def XYZCoord(self, eta): 
        """
        Parameters
        ----------
        eta

        Returns
        -------
        XYZ
        """
        # if('1' in self.regType):
            # xCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,0],eta,method='linear')
            # yCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,1],eta,method='linear')
            # zCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,2],eta,method='linear')
            # return np.transpose(np.array([xCrd,yCrd,zCrd]))
        # elif('2' in self.regType or '3' in self.regType):
            # xCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,0],eta,method='cubic')
            # yCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,1],eta,method='cubic')
            # zCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,2],eta,method='cubic')
            # return np.transpose(np.array([xCrd,yCrd,zCrd]))
        numPts = len(eta)
        if('quad1' == self.regType):
            #Nvec = np.zeros((1,4))
            Nmat = np.zeros((numPts,4))
            for i in range(0,numPts):
                Nmat[i,0] = 0.25*(eta[i,0] - 1.0)*(eta[i,1] - 1.0)
                Nmat[i,1] = -0.25*(eta[i,0] + 1.0)*(eta[i,1] - 1.0)
                Nmat[i,2] = 0.25*(eta[0] + 1.0)*(eta[i,1] + 1.0)
                Nmat[i,3] = -0.25*(eta[0] - 1.0)*(eta[i,1] + 1.0)
            XYZ = np.matmul(Nmat,self.keyPts)
        elif('quad2' == self.regType):
            r1 = - 1
            r2 = 0
            r3 = 1
            Nmat = np.zeros((numPts,9))
            for i in range(0,numPts):
                Nmat[i,0] = 0.25*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,1] - r2)*(eta[i,1] - r3)
                Nmat[i,1] = 0.25*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,1] - r2)*(eta[i,1] - r3)
                Nmat[i,2] = 0.25*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,1] - r1)*(eta[i,1] - r2)
                Nmat[i,3] = 0.25*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,1] - r1)*(eta[i,1] - r2)
                Nmat[i,4] = - 0.5*(eta[i,0] - r1)*(eta[i,0] - r3)*(eta[i,1] - r2)*(eta[i,1] - r3)
                Nmat[i,5] = - 0.5*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,1] - r1)*(eta[i,1] - r3)
                Nmat[i,6] = - 0.5*(eta[i,0] - r1)*(eta[i,0] - r3)*(eta[i,1] - r1)*(eta[i,1] - r2)
                Nmat[i,7] = - 0.5*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,1] - r1)*(eta[i,1] - r3)
                Nmat[i,8] = (eta[i,0] - r1)*(eta[i,0] - r3)*(eta[i,1] - r1)*(eta[i,1] - r3)
            XYZ = np.matmul(Nmat,self.keyPts)
        elif('quad3' == self.regType):
            r1 = - 1
            r2 = - 0.333333333333333
            r3 = 0.333333333333333
            r4 = 1
            coef = np.array([0.31640625,- 0.31640625,0.31640625,- 0.31640625,- 0.94921875,0.94921875,0.94921875,- 0.94921875,- 0.94921875,0.94921875,0.94921875,- 0.94921875,2.84765625,- 2.84765625,2.84765625,- 2.84765625])
            Nmat = np.zeros((numPts,16))
            for i in range(0,numPts):
                Nmat[i,0] = coef[0]*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,0] - r4)*(eta[i,1] - r2)*(eta[i,1] - r3)*(eta[i,1] - r4)
                Nmat[i,1] = coef[1]*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,1] - r2)*(eta[i,1] - r3)*(eta[i,1] - r4)
                Nmat[i,2] = coef[2]*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,1] - r1)*(eta[i,1] - r2)*(eta[i,1] - r3)
                Nmat[i,3] = coef[3]*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,0] - r4)*(eta[i,1] - r1)*(eta[i,1] - r2)*(eta[i,1] - r3)
                Nmat[i,4] = coef[4]*(eta[i,0] - r1)*(eta[i,0] - r3)*(eta[i,0] - r4)*(eta[i,1] - r2)*(eta[i,1] - r3)*(eta[i,1] - r4)
                Nmat[i,5] = coef[5]*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,0] - r4)*(eta[i,1] - r2)*(eta[i,1] - r3)*(eta[i,1] - r4)
                Nmat[i,6] = coef[6]*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,1] - r1)*(eta[i,1] - r3)*(eta[i,1] - r4)
                Nmat[i,7] = coef[7]*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,1] - r1)*(eta[i,1] - r2)*(eta[i,1] - r4)
                Nmat[i,8] = coef[8]*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,0] - r4)*(eta[i,1] - r1)*(eta[i,1] - r2)*(eta[i,1] - r3)
                Nmat[i,9] = coef[9]*(eta[i,0] - r1)*(eta[i,0] - r3)*(eta[i,0] - r4)*(eta[i,1] - r1)*(eta[i,1] - r2)*(eta[i,1] - r3)
                Nmat[i,10] = coef[10]*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,0] - r4)*(eta[i,1] - r1)*(eta[i,1] - r2)*(eta[i,1] - r4)
                Nmat[i,11] = coef[11]*(eta[i,0] - r2)*(eta[i,0] - r3)*(eta[i,0] - r4)*(eta[i,1] - r1)*(eta[i,1] - r3)*(eta[i,1] - r4)
                Nmat[i,12] = coef[12]*(eta[i,0] - r1)*(eta[i,0] - r3)*(eta[i,0] - r4)*(eta[i,1] - r1)*(eta[i,1] - r3)*(eta[i,1] - r4)
                Nmat[i,13] = coef[13]*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,0] - r4)*(eta[i,1] - r1)*(eta[i,1] - r3)*(eta[i,1] - r4)
                Nmat[i,14] = coef[14]*(eta[i,0] - r1)*(eta[i,0] - r2)*(eta[i,0] - r4)*(eta[i,1] - r1)*(eta[i,1] - r2)*(eta[i,1] - r4)
                Nmat[i,15] = coef[15]*(eta[i,0] - r1)*(eta[i,0] - r3)*(eta[i,0] - r4)*(eta[i,1] - r1)*(eta[i,1] - r2)*(eta[i,1] - r4)
            XYZ = np.matmul(Nmat,self.keyPts)
        elif('tri1' == self.regType):
            Nmat = np.zeros((numPts,3))
            for i in range(0,numPts):
                Nmat[i,0] = 1 - eta[i,0] - eta[i,1]
                Nmat[i,1] = eta[i,0]
                Nmat[i,2] = eta[i,1]
            XYZ = np.matmul(Nmat,self.keyPts)
        elif('tri2' == self.regType):
            Nmat = np.zeros((numPts,6))
            for i in range(0,numPts):
                Nmat[i,0] = 2*(eta[i,0] + eta[i,1] - 1)*(eta[i,0] + eta[i,1] - 0.5)
                Nmat[i,1] = 2*eta[i,0]*(eta[i,0] - 0.5)
                Nmat[i,2] = 2*eta[i,1]*(eta[i,1] - 0.5)
                Nmat[i,3] = - 4*eta[i,0]*(eta[i,0] + eta[i,1] - 1)
                Nmat[i,4] = 4*eta[i,0]*eta[i,1]
                Nmat[i,5] = - 4*eta[i,1]*(eta[i,0] + eta[i,1] - 1)
            XYZ = np.matmul(Nmat,self.keyPts)
        elif 'tri3' == self.type:
            r2 = 1 / 3
            r3 = 2 / 3
            coef = np.array([- 4.5,4.5,4.5,13.5,- 13.5,13.5,13.5,- 13.5,13.5,- 27])
            Nmat = np.zeros((numPts,10))
            for i in range(0,numPts):
                Nmat[i,0] = coef[0]*(eta[i,0] + eta[i,1] - r2)*(eta[i,0] + eta[i,1] - r3)*(eta[i,0] + eta[i,1] - 1)
                Nmat[i,1] = coef[1]*eta[i,0]*(eta[i,0] - r2)*(eta[i,0] - r3)
                Nmat[i,2] = coef[2]*eta[i,1]*(eta[i,1] - r2)*(eta[i,1] - r3)
                Nmat[i,3] = coef[3]*eta[i,0]*(eta[i,0] + eta[i,1] - r3)*(eta[i,0] + eta[i,1] - 1)
                Nmat[i,4] = coef[4]*eta[i,0]*(eta[i,0] - r2)*(eta[i,0] + eta[i,1] - 1)
                Nmat[i,5] = coef[5]*eta[i,0]*eta[i,1]*(eta[i,0] - r2)
                Nmat[i,6] = coef[6]*eta[i,0]*eta[i,1]*(eta[i,1] - r2)
                Nmat[i,7] = coef[7]*eta[i,1]*(eta[i,1] - r2)*(eta[i,0] + eta[i,1] - 1)
                Nmat[i,8] = coef[8]*eta[i,1]*(eta[i,0] + eta[i,1] - r3)*(eta[i,0] + eta[i,1] - 1)
                Nmat[i,9] = coef[9]*eta[i,0]*eta[i,1]*(eta[i,0] + eta[i,1] - 1)
            XYZ = np.matmul(Nmat,self.keyPts)
        elif(self.regType == 'sphere'):
            vec = self.keyPts[1,:] - self.keyPts[0,:] ## local x-direction
            outerRad = np.linalg.norm(vec)
            XYZ = np.zeros((len(eta),3))
            ndi = 0
            for nd in eta:
                phiComp =  np.linalg.norm(nd)
                if (nd[0] > 1.0e-12):
                    theta = np.arctan(nd[1]/nd[0])
                else:
                    if (nd[0] < 1.0e-12):
                        theta = np.pi + np.arctan(nd[1] / nd[0])
                    else:
                        theta = np.arctan(nd[1]/1.0e-12)
                phi = 0.5*np.pi - phiComp
                xloc = outerRad*np.cos(theta)*np.cos(phi)
                yloc = outerRad*np.sin(theta)*np.cos(phi)
                zloc = outerRad*np.sin(phi)
                XYZLoc = np.array([xloc,yloc,zloc])
                a1 = (1/outerRad)*vec
                vec2 = self.keyPts[2,:] - self.keyPts[1,:]
                vec3 = np.array([
                    (vec[1]*vec2[2] - vec[2]*vec2[1]),
                    (vec[2]*vec2[0] - vec[0]*vec2[2]),
                    (vec[0]*vec2[1] - vec[1]*vec2[0])
                    ])
                mag = np.sqrt(vec3*vec3.T)
                a3 = (1/mag)*vec3
                a2 = np.array([
                    (a3[1]*a1[2] - a3[2]*a1[1]),
                    (a3[2]*a1[0] - a3[0]*a1[2]),
                    (a3[0]*a1[1] - a3[1]*a1[0])
                    ])
                alpha = np.array([[a1],[a2],[a3]])
                XYZ[ndi] = np.matmul(XYZLoc,alpha) + self.keyPts[0,:]
                ndi = ndi + 1
        
        return XYZ