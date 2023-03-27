import numpy as np

class SpatialGridList3D():

    def __init__(self, minimumX, maximumX, minimumY, maximumY, minimumZ, maximumZ, xGridSize, yGridSize, zGridSize):
        self.xMin = minimumX
        self.yMin = minimumY
        self.zMin = minimumZ
        self.xGSz = xGridSize
        self.yGSz = yGridSize
        self.zGSz = zGridSize
        xLen = maximumX - minimumX
        yLen = maximumY - minimumY
        zLen = maximumZ - minimumZ
        self.xRows = int(np.ceil(xLen/xGridSize))
        self.yRows = int(np.ceil(yLen/yGridSize))
        self.zRows = int(np.ceil(zLen/zGridSize))
        self.fullList = list()
        for i in range(0,self.xRows):
            xList = list()
            for j in range(0,self.yRows):
                yList = list()
                for k in range(0,self.zRows):
                    yList.append(list())
                xList.append(yList)
            self.fullList.append(xList)
            
    def addEntry(self, val, coord):
        xRow = int(np.floor((coord[0] - self.xMin)/self.xGSz))
        yRow = int(np.floor((coord[1] - self.yMin)/self.yGSz))
        zRow = int(np.floor((coord[2] - self.zMin)/self.zGSz))
        self.fullList[xRow][yRow][zRow].append(val)
        
    def findInXYZMargin(self,point,Xmargin,Ymargin,Zmargin):
        if(Xmargin == -1):
            iMax = self.xRows
            iMin = 0
        else:
            iMax = int(np.ceil((point[0] + Xmargin - self.xMin)/self.xGSz))
            if(iMax > self.xRows):
                iMax = self.xRows
            iMin = int(np.floor((point[0] - Xmargin - self.xMin)/self.xGSz))
            if(iMin < 0):
                iMin = 0
                
        if(Ymargin == -1):
            jMax = self.yRows
            jMin = 0
        else:
            jMax = int(np.ceil((point[1] + Ymargin - self.yMin)/self.yGSz))
            if(jMax > self.yRows):
                jMax = self.yRows
            jMin = int(np.floor((point[1] - Ymargin - self.yMin)/self.yGSz))
            if(jMin < 0):
                jMin = 0
        
        if(Zmargin == -1):
            kMax = self.zRows
            kMin = 0
        else:
            kMax = int(np.ceil((point[2] + Zmargin - self.zMin)/self.zGSz))
            if(kMax > self.zRows):
                kMax = self.zRows
            kMin = int(np.floor((point[2] - Zmargin - self.zMin)/self.zGSz))
            if(kMin < 0):
                kMin = 0
        
        labelList = list()
        for i in range(iMin,iMax):
            for j in range(jMin,jMax):
                for k in range(kMin,kMax):
                    labelList.extend(self.fullList[i][j][k])
        
        return labelList
        
    def findInRadius(self, point, radius):
        labelList = self.findInXYZMargin(point,radius,radius,radius)
        return labelList