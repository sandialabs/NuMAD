classdef elementSet
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        plygroups
        elementList
    end
    
    methods
        function obj = elementSet(setName,setPlyGroups,setElList)
            obj.name = setName;
            obj.plygroups = setPlyGroups;
            obj.elementList = setElList;
        end
    end
end

