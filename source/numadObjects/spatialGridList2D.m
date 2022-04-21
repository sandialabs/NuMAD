classdef spatialGridList2D
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        firstEnt = [];
        label = [];
        nextEnt = [];
        xMin = 0;
        yMin = 0;
        xGSz = 1;
        yGSz = 1;
        xRows = 1;
        yRows = 1;
    end
    
    methods
        function obj = spatialGridList2D(minimumX,maximumX,minimumY,maximumY,xGridSize,yGridSize)
            obj.xMin = minimumX;
            obj.yMin = minimumY;
            obj.xGSz = xGridSize;
            obj.yGSz = yGridSize;
            XLen = maximumX - minimumX;
            YLen = maximumY - minimumY;
            obj.xRows = ceil(XLen/xGridSize);
            obj.yRows = ceil(YLen/yGridSize);
            obj.firstEnt = zeros(obj.xRows,obj.yRows);
            obj.label = [];
            obj.nextEnt = [];
        end
        
        function obj = addEntry(obj,val,coord)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xRow = ceil((coord(1) - obj.xMin)/obj.xGSz);
            yRow = ceil((coord(2) - obj.yMin)/obj.xGSz);
            if(obj.firstEnt(xRow,yRow) == 0)
                obj.label = [obj.label;val];
                obj.nextEnt = [obj.nextEnt;0];
                obj.firstEnt(xRow,yRow) = length(obj.label);
            else
                i = obj.firstEnt(xRow,yRow);
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
            labelList = obj.findInXYMargin(point,radius,radius);
        end
        
        function labelList = findInXYMargin(obj,point,Xmargin,Ymargin)
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
            
            labelList = [];
            for i = iMin:iMax
                for j = jMin:jMax
                    if(obj.firstEnt(i,j) ~= 0)
                        k = obj.firstEnt(i,j);
                        endReached = 0;
                        while(endReached == 0)
                            labelList = [labelList,obj.label(k)];
                            if(obj.nextEnt(k) == 0)
                                endReached = 1;
                            else
                                k = obj.nextEnt(k);
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

