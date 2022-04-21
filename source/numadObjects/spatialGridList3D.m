classdef spatialGridList3D
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        firstEnt = [];
        label = [];
        nextEnt = [];
        xMin = 0;
        yMin = 0;
        zMin = 0;
        xGSz = 1;
        yGSz = 1;
        zGSz = 1;
        xRows = 1;
        yRows = 1;
        zRows = 1;
    end
    
    methods
        function obj = spatialGridList3D(minimumX,maximumX,minimumY,maximumY,minimumZ,maximumZ,xGridSize,yGridSize,zGridSize)
            obj.xMin = minimumX;
            obj.yMin = minimumY;
            obj.zMin = minimumZ;
            obj.xGSz = xGridSize;
            obj.yGSz = yGridSize;
            obj.zGSz = zGridSize;
            XLen = maximumX - minimumX;
            YLen = maximumY - minimumY;
            ZLen = maximumZ - minimumZ;
            obj.xRows = ceil(XLen/xGridSize);
            obj.yRows = ceil(YLen/yGridSize);
            obj.zRows = ceil(ZLen/zGridSize);
            obj.firstEnt = zeros(obj.xRows,obj.yRows,obj.zRows);
            obj.label = [];
            obj.nextEnt = [];
        end
        
        function obj = addEntry(obj,val,coord)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xRow = ceil((coord(1) - obj.xMin)/obj.xGSz);
            yRow = ceil((coord(2) - obj.yMin)/obj.yGSz);
            zRow = ceil((coord(3) - obj.zMin)/obj.zGSz);
            if(obj.firstEnt(xRow,yRow,zRow) == 0)
                obj.label = [obj.label;val];
                obj.nextEnt = [obj.nextEnt;0];
                obj.firstEnt(xRow,yRow,zRow) = length(obj.label);
            else
                i = obj.firstEnt(xRow,yRow,zRow);
                inserted = 0;
                while(inserted == 0)
                    if(obj.nextEnt(i) == 0)
                        obj.label = [obj.label;val];
                        obj.nextEnt = [obj.nextEnt;0];
                        obj.nextEnt(i) = length(obj.label);
                        inserted = 1;
                    else
                        i = obj.nextEnt(i);
                    end
                end
            end
        end
        
        function [labelList] = findInRadius(obj,point,radius)
            labelList = obj.findInXYZMargin(point,radius,radius,radius);
        end
        
        function labelList = findInXYZMargin(obj,point,Xmargin,Ymargin,Zmargin)
            if(Xmargin == -1)
                iMax = obj.xRows;
                iMin = 1;
            else
                iMax = ceil((point(1) + Xmargin - obj.xMin)/obj.xGSz);
                iMax = min([iMax,obj.xRows],[],'all');
                iMin = ceil((point(1) - Xmargin - obj.xMin)/obj.xGSz);
                iMin = max([iMin,1],[],'all');
            end
            if(Ymargin == -1)
                jMax = obj.yRows;
                jMin = 1;
            else
                jMax = ceil((point(2) + Ymargin - obj.yMin)/obj.yGSz);
                jMax = min([jMax,obj.yRows],[],'all');
                jMin = ceil((point(2) - Ymargin - obj.yMin)/obj.yGSz);
                jMin = max([jMin,1],[],'all');
            end
            if(Zmargin == -1)
                kMax = obj.zRows;
                kMin = 1;
            else
                kMax = ceil((point(3) + Zmargin - obj.zMin)/obj.zGSz);
                kMax = min([kMax,obj.zRows],[],'all');
                kMin = ceil((point(3) - Zmargin - obj.zMin)/obj.zGSz);
                kMin = max([kMin,1],[],'all');
            end
            
            labelList = [];
            for i = iMin:iMax
                for j = jMin:jMax
                    for k = kMin:kMax
                        if(obj.firstEnt(i,j,k) ~= 0)
                            m = obj.firstEnt(i,j,k);
                            endReached = 0;
                            while(endReached == 0)
                                labelList = [labelList,obj.label(m)];
                                if(obj.nextEnt(m) == 0)
                                    endReached = 1;
                                else
                                    m = obj.nextEnt(m);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function obj = reOrderLabels(obj,labelOrder)
            currentLabel = zeros(max(labelOrder,[],'all'),1);
            for i = 1:length(labelOrder)
                currentLabel(labelOrder(i)) = i;
            end
            for i = 1:length(obj.label)
                orig = obj.label(i);
                if(orig ~= 0)
                    obj.label(i) = currentLabel(orig);
                end
            end
        end
    end
end

