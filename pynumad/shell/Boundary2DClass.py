import numpy as np
import pynumad.shell.MeshTools as mt
from pynumad.shell.Segment2DClass import *

class Boundary2D():

    def __init__(self,segList=[]):
        self.segList = list()
        self.segList.extend(segList)
        
    def addSegment(self,segType,keyPts,numEls):
        self.segList.append(Segment2D(segType,keyPts,numEls))
        
    def getBoundaryMesh(self):
        allNds = list()
        allEds = list()
        totNds = 0
        for seg in self.segList:
            segMesh = seg.getNodesEdges()
            allNds.extend(segMesh['nodes'])
            allEds.extend(segMesh['edges'] + totNds)
            totNds = len(allNds)
        allNds = np.array(allNds)
        allEds = np.array(allEds)
        
        meshData = dict()
        meshData['nodes'] = allNds
        meshData['elements'] = allEds
        
        output = mt.mergeDuplicateNodes(meshData)
        
        return output