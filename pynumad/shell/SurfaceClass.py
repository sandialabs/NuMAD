import numpy as np
import pynumad.shell.MeshTools as mt
from pynumad.shell.ShellRegionClass import *

class Surface():

    def __init__(self,regionList=[],regionNames=[],meshList=[],meshNames=[]):
        self.shellRegions = list()
        self.shellRegions.extend(regionList)
        self.regionNames = list()
        self.regionNames.extend(regionNames)
        self.meshes = list()
        self.meshes.extend(meshList) 
        self.meshNames = list()
        self.meshNames.extend(meshNames)
        
    def addShellRegion(self,regType,keyPts,numEls,name=None,natSpaceCrds=[],elType='quad',meshMethod='free'):
        self.shellRegions.append(ShellRegion(regType,keyPts,numEls,natSpaceCrds,elType,meshMethod))
        if(name == None):
            numReg = len(self.shellRegions)
            regName = 'Sub-Region_' + str(numReg)
            self.regionNames.append(regName)
        else:
            self.regionNames.append(name)
        
    def addMesh(self,meshData,name=None):
        self.meshes.append(meshData)
        if(name == None):
            numMsh = len(self.meshes)
            meshName = 'Sub-Mesh_' + str(numMsh)
            self.meshNames.append(meshName)
        else:
            self.meshNames.append(name)
        
    def getSurfaceMesh(self):
        allNds = list()
        allEls = list()
        elSetList = list()
        numNds = 0
        numEls = 0
        regi = 0
        for reg in self.shellRegions:
            regMesh = reg.createShellMesh()
            allNds.extend(regMesh['nodes'])
            setList = list()
            eli = 0
            for el in regMesh['elements']:
                for i in range(0,4):
                    if(el[i] != -1):
                        el[i] = el[i] + numNds
                allEls.append(el)
                setList.append((eli + numEls))
                eli = eli + 1
            thisSet = dict()
            thisSet['name'] = self.regionNames[regi]
            thisSet['labels'] = setList
            elSetList.append(thisSet)
            numNds = len(allNds)
            numEls = len(allEls)
            regi = regi + 1
        mshi = 0
        for msh in self.meshes:
            allNds.extend(msh['nodes'])
            setList = list()
            eli = 0
            for el in msh['elements']:
                for i in range(0,4):
                    if(el[i] != -1):
                        el[i] = el[i] + numNds
                allEls.append(el)
                setList.append((eli + numEls))
                eli = eli + 1
            thisSet = dict()
            thisSet['name'] = self.meshNames[mshi]
            thisSet['labels'] = setList
            elSetList.append(thisSet)
            numNds = len(allNds)
            numEls = len(allEls)
            mshi = mshi + 1
        mData = dict()
        mData['nodes'] = np.array(allNds)
        mData['elements'] = np.array(allEls)
        mData = mt.mergeDuplicateNodes(mData)
        mData['sets'] = dict()
        mData['sets']['element'] = elSetList
        return mData