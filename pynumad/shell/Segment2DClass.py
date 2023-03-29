import numpy as np
from scipy import interpolate

class Segment2D():

    def __init__(self,segType,keyPts,numEls):
        self.segType = segType ## line, curve, arc
        self.keyPts = keyPts
        self.numEls = numEls
        
    def getNodesEdges(self):
        nNds = self.numEls+1
        if(self.segType == 'line'):
            pt1 = np.array(self.keyPts[0])
            pt2 = np.array(self.keyPts[1])
            proj = pt2 - pt1
            steps = (1.0/self.numEls)*np.array(range(0,nNds))
            nds = list()
            for st in steps:
                nd = pt1 + st*proj
                nds.append(nd)
            nodes = np.array(nds)
            eN1 = np.array(range(0,self.numEls),dtype=int)
            eN2 = np.array(range(1,nNds),dtype=int)
            edges = np.transpose(np.array([eN1,eN2]))
            output = dict()
            output['nodes'] = nodes
            output['edges'] = edges
            return output
        elif(self.segType == 'curve'):
            kPTp = np.transpose(np.array(self.keyPts))
            numKp = len(self.keyPts)
            pKp = (1.0/(numKp-1))*np.array(range(0,numKp))
            pNds = (1.0/self.numEls)*np.array(range(0,nNds))
            if(numKp == 2):
                order = 'linear'
            elif(numKp == 3):
                order = 'quadratic'
            else:
                order = 'cubic'
            iFun = interpolate.interp1d(pKp,kPTp[0],order, axis=0,bounds_error=False,fill_value='extrapolate')
            xNds = iFun(pNds)
            iFun = interpolate.interp1d(pKp,kPTp[1],order, axis=0,bounds_error=False,fill_value='extrapolate')
            yNds = iFun(pNds)
            eN1 = np.array(range(0,self.numEls),dtype=int)
            eN2 = np.array(range(1,nNds),dtype=int)
            output = dict()
            output['nodes'] = np.transpose(np.array([xNds,yNds]))
            output['edges'] = np.transpose(np.array([eN1,eN2]))
            return output
        elif(self.segType == 'arc'):
            kPar = np.array(self.keyPts)
            dist13 = np.linalg.norm(kPar[0] - kPar[2])
            if(dist13 < 1.0e-12):
                center = 0.5*(kPar[0] + kPar[1])
                rad = np.linalg.norm(kPar[0] - center)
            else:
                center = 0.3333*(kPar[0] + kPar[1] + kPar[2])
                ndc = 1.0
                i = 0
                Rvec = np.zeros(2)
                dRdC = np.zeros((2,2))
                while(ndc > 1.0e-12 and i < 50):
                    v1 = kPar[0] - center
                    v2 = kPar[1] - center
                    v3 = kPar[2] - center
                    Rvec[0] = np.dot(v1,v1) - np.dot(v2,v2)
                    Rvec[1] = np.dot(v1,v1) - np.dot(v3,v3)
                    dRdC[0] = 2.0*(v2 - v1)
                    dRdC[1] = 2.0*(v3 - v1)
                    dc = np.linalg.solve(dRdC,-Rvec)
                    center = center + dc
                    ndc = np.linalg.norm(dc)
                    i = i + 1
                rad = np.linalg.norm(kPar[0] - center)
            theta = np.zeros(3)
            for i in range(0,3):
                xRel = kPar[i,0] - center[0]
                yRel = kPar[i,1] - center[1]
                if(abs(xRel) < 1.0e-12):
                    xRel = 1.0e-12
                if(xRel > 0.0):
                    theta[i] = np.arctan(yRel/xRel)
                else:
                    theta[i] = np.arctan(yRel/xRel) + np.pi
            if(theta[1] > theta[0] and theta[1] < theta[2]):
                thetaNds = np.linspace(theta[0],theta[2],nNds)
            else:
                t0adj = theta[0] + 2.0*np.pi
                thetaNds = np.linspace(t0adj,theta[2],nNds)
            nodes = np.zeros((nNds,2))
            for thi in range(0,nNds):
                nodes[thi,0] = rad*np.cos(thetaNds[thi]) + center[0]
                nodes[thi,1] = rad*np.sin(thetaNds[thi]) + center[1]
            eN1 = np.array(range(0,self.numEls),dtype=int)
            eN2 = np.array(range(1,nNds),dtype=int)
            output = dict()
            output['nodes'] = nodes
            output['edges'] = np.transpose(np.array([eN1,eN2]))
            return output