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

                XYZ = self.XYZCoord(mData['nodes'])
                
                mData['nodes'] = XYZ
                mData['elements'] = elLst
                return mData

            else:
                raise Exception('Only quadrilateral shell regions can use the structured meshing option')

        else:
            bndData = self.initialBoundary()
            mesh = Mesh2D(bndData['nodes'],bndData['elements'])
            # print('boundaryElements')
            # print(bndData['elements'])
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
        if('1' in self.regType):
            xCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,0],eta,method='linear')
            yCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,1],eta,method='linear')
            zCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,2],eta,method='linear')
            return np.transpose(np.array([xCrd,yCrd,zCrd]))
        elif('2' in self.regType or '3' in self.regType):
            xCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,0],eta,method='cubic')
            yCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,1],eta,method='cubic')
            zCrd = interpolate.griddata(self.natSpaceCrd,self.keyPts[:,2],eta,method='cubic')
            return np.transpose(np.array([xCrd,yCrd,zCrd]))
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
                xloc = outerRad * np.cos(theta) * np.cos(phi)
                yloc = outerRad * np.sin(theta) * np.cos(phi)
                zloc = outerRad * np.sin(phi)
                XYZLoc = np.array([xloc,yloc,zloc])
                a1 = (1/outerRad)*vec
                vec2 = self.keyPts[2,:] - self.keyPts[1,:]
                vec3 = np.array([
                    (vec[1] * vec2[2] - vec[2] * vec2[1]),
                    (vec[2] * vec2[0] - vec[0] * vec2[2]),
                    (vec[0] * vec2[1] - vec[1] * vec2[0])
                    ])
                mag = np.sqrt(vec3 * vec3.T)
                a3 = (1/mag)*vec3
                a2 = np.array([
                    (a3[1] * a1[2] - a3[2] * a1[1]),
                    (a3[2] * a1[0] - a3[0] * a1[2]),
                    (a3[0] * a1[1] - a3[1] * a1[0])
                    ])
                alpha = np.array([[a1],[a2],[a3]])
                XYZ[ndi] = np.matmul(XYZLoc,alpha) + self.keyPts[0,:]
                ndi = ndi + 1
            return XYZ