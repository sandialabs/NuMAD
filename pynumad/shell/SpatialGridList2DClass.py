import numpy as np

class SpatialGridList2D():

    def __init__(self, minimumX, maximumX, minimumY, maximumY, xGridSize, yGridSize):
        self.xMin = minimumX
        self.yMin = minimumY
        self.xGSz = xGridSize
        self.yGSz = yGridSize
        xLen = maximumX - minimumX
        yLen = maximumY - minimumY
        self.xRows = int(np.ceil(xLen/xGridSize))
        self.yRows = int(np.ceil(yLen/yGridSize))
        self.fullList = list()
        for i in range(0,self.xRows):
            xList = list()
            for j in range(0,self.yRows):
                xList.append(list())
            self.fullList.append(xList)
    
    def addEntry(self, val, coord):
        xRow = int(np.floor((coord[0] - self.xMin)/self.xGSz))
        yRow = int(np.floor((coord[1] - self.yMin)/self.yGSz))
        self.fullList[xRow][yRow].append(val)
        
    def findInXYMargin(self,point,Xmargin,Ymargin):
        if(Xmargin == -1):
            iMax = self.xRows
            iMin = 0
        else:
            # outStr = 'point[0] ' + str(point[0]) + ' Xmargin ' + str(Xmargin) + ' xMin ' + str(self.xMin) + ' xGSz ' + str(self.xGSz)
            # print(outStr)
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
        
        labelList = list()
        for i in range(iMin,iMax):
            for j in range(jMin,jMax):
                labelList.extend(self.fullList[i][j])
        
        return labelList
        
    def findInRadius(self, point, radius):
        labelList = self.findInXYMargin(point,radius,radius)
        return labelList