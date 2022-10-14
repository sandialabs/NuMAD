classdef NuMesh3D
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nodes = [];
        elements = [];
        boundaryFaces = [];
    end
    
    methods
        function obj = NuMesh3D(boundaryNodes,boundaryFaces)
            obj.nodes = boundaryNodes;
            obj.boundaryFaces = boundaryFaces;
        end
        
        function [nodes,elements,obj] = createSweptMesh(obj,sweepMethod,direction,sweepDistance,sweepElements,followNormal,destNodes)
            if(contains(sweepMethod,'to_dest_nodes'))
                totalSwpEls = sum(sweepElements,'all');
                guideNodes = [];
                for j = 1:size(obj.nodes,1)
                    guideNodes = [guideNodes;obj.nodes(j,:)'];
                end
                destLen = size(destNodes,2);
                for i = 1:destLen
                    Dnds = destNodes{i};
                    col = [];
                    for j = 1:size(obj.nodes,1)
                        col = [col;Dnds(j,:)'];
                    end
                    guideNodes = [guideNodes,col];
                end
                allNodes = [];
                guideParVal = 0;
                for i = 1:length(sweepElements)
                    nextEnt = guideParVal(end) + sweepElements(i);
                    guideParVal = [guideParVal,nextEnt];
                end
                guideParVal = (1/nextEnt)*guideParVal;
%                 guideParVal = linspace(0,1,destLen+1);
                allParVal = linspace(0,1,totalSwpEls+1);
                for i = 1:size(guideNodes,1)
                    row = interp1(guideParVal,guideNodes(i,:),allParVal,'pchip');
                    allNodes = [allNodes;row];
                end
                nodes = [];
                for i = 1:totalSwpEls+1
                    for j = 1:3:size(allNodes,1)
                        nodes = [nodes;allNodes(j:j+2,i)'];
                    end
                end
                numLayerNds = size(obj.nodes,1);
                numLayerEls = size(obj.boundaryFaces,1);
                elements = [];
                four1s = ones(1,4);
                for i = 1:totalSwpEls
                    for j = 1:numLayerEls
                        face = obj.boundaryFaces(j,:);
                        lowFace = face + numLayerNds*(i-1)*four1s;
                        highFace = face + numLayerNds*i*four1s;
                        if(face(4) == 0)
                            newEl = [lowFace(1:3),highFace(1:3),0,0];
                            n1 = newEl(1);
                            n2 = newEl(2);
                            n3 = newEl(3);
                            n4 = newEl(4);
                            v1 = nodes(n2,:) - nodes(n1,:);
                            v2 = nodes(n3,:) - nodes(n1,:);
                            v3 = nodes(n4,:) - nodes(n1,:);
                            vMat = [v1;v2;v3];
                            detV = det(vMat);
                            if(detV < 0)
                                swap = newEl(2);
                                newEl(2) = newEl(3);
                                newEl(3) = swap;
                                swap = newEl(5);
                                newEl(5) = newEl(6);
                                newEl(6) = swap;
                            end
                        else
                            newEl = [lowFace,highFace];
                            n1 = newEl(1);
                            n2 = newEl(2);
                            n3 = newEl(3);
                            n5 = newEl(5);
                            v1 = nodes(n2,:) - nodes(n1,:);
                            v2 = nodes(n3,:) - nodes(n1,:);
                            v3 = nodes(n5,:) - nodes(n1,:);
                            vMat = [v1;v2;v3];
                            detV = det(vMat);
                            if(detV < 0)
                                swap = newEl(2);
                                newEl(2) = newEl(4);
                                newEl(4) = swap;
                                swap = newEl(6);
                                newEl(6) = newEl(8);
                                newEl(8) = swap;
                            end
                        end
                        elements = [elements;newEl];
                    end
                end
                obj.nodes = nodes;
                obj.elements = elements;
            end
        end
        
    end
end

