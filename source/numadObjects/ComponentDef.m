%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Part of the SNL NuMAD Toolbox                    
%  Developed by Sandia National Laboratories Wind Energy Technologies 
%              See license.txt for disclaimer information             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef ComponentDef < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ``ComponentDef``  A class definition for blade components.
%
% Examples: 
% 
%	``comp_obj = ComponentDef();``
% 
%	``comp_obj = ComponentDef(comp_struct);``
%
% See also ``xlsBlade``, ``BladeDef``, ``BladeDef.addComponent``
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        group                   % Integar: 0 = blade, 1 = first shear web, 2 = second shear web, etc.
        name = ''               % String: Name, such as 'spar'
        materialid              % String: Material id number from blade.materials
        fabricangle             % Float: Fiber angle
        hpextents               % String Array: Array of keypoints such as {'b','c'}
        lpextents               % String Array: Array of keypoints such as {'b','c'}
        cp                      % control points defining layer distribution
        imethod = 'pchip'       % String: imethod
        pinnedends = false;     %
    end
    properties (Hidden)
        hCtrl                   % ?
        hLine                   % ?
    end
    
    methods 
        function obj = ComponentDef(comp)
            if nargin > 0
                obj.group        = comp.group;
                obj.name         = comp.name;
                obj.materialid   = comp.materialid;
                obj.fabricangle  = comp.fabricangle;
                obj.hpextents    = comp.hpextents;
                obj.lpextents    = comp.lpextents;
                obj.cp           = comp.cp;
                obj.imethod      = comp.imethod;
            end
        end
        
        function [cpx,cpy] = getcp(obj)
            if obj.pinnedends
                if any(obj.cp(:,1)<0) || any(obj.cp(:,1)>1)
                    error('ComponentDef: first coordinate of control points must be in range [0,1] when using "pinned" ends');
                end
                cpx = [-0.01; obj.cp(:,1); 1.01];
                cpy = [0; obj.cp(:,2); 0];
            else
                cpx = obj.cp(:,1);
                cpy = obj.cp(:,2);
            end
        end
        
        function nLayers = getNumLayers(obj,span)
            [cpx,cpy] = getcp(obj);
            try
                nLayers = interp1(cpx,cpy,span,obj.imethod,0);
            catch
                nLayers = interp1(cpx,cpy,span,'linear',0);
            end
        end
        
        function plotcp(obj)
            [cpx,cpy] = getcp(obj);
            obj.hCtrl = line(cpx,cpy,'Marker','s','LineStyle','none');
%             x = linspace(min(cpx),max(cpx),100);
            x = linspace(0,1,100);
            y = round(interp1(cpx,cpy,x,'pchip',0));
            obj.hLine = line(x,y,'LineStyle','--');
            set(obj.hCtrl,'ButtonDownFcn',@bdf_movepts,'UserData',obj);
            set(obj.hLine,'HitTest','off');
            title(obj.name);
        end
    end
end

function bdf_movepts(cbo,~)
    click = get(gca,'CurrentPoint');  %get pointer position within current axes
    obj = get(cbo,'UserData');
    pt = click(1,1:2);     %get x-y values from pointer position
    [cpx,cpy] = getcp(obj);
    dist = hypot(cpx-pt(1),cpy-pt(2));
    [m, i] = min(dist);  % find the closest point
    if obj.pinnedends && any(i==[1,numel(cpx)])
        return;  % do not allow moving the "pinned" end points
    end
    ax = axis;
    if m < 5e-2*max((ax(2)-ax(1)),(ax(4)-ax(3)))
        rbbox;
        release = get(gca,'CurrentPoint');
        pt = release(1,1:2);
        if obj.pinnedends
            pt(1) = min(max(pt(1),0),1);  % constrain x to [0,1]
        end
        pt(2) = max(pt(2),0);  % constrain y >= 0
        if i>1, pt(1) = max(pt(1),cpx(i-1)+eps(cpx(i-1))); end  % don't allow movement left of another cp
        if i<numel(cpx), pt(1) = min(pt(1),cpx(i+1)-eps(cpx(i+1))); end % don't allow movement right of another cp
        cpx(i) = pt(1);
        cpy(i) = pt(2);
        set(obj.hCtrl,'XData',cpx,'YData',cpy);
        x = get(obj.hLine,'XData');
        y = round(interp1(cpx,cpy,x,'pchip',0));
        set(obj.hLine,'XData',x,'YData',y);
        if obj.pinnedends
            obj.cp(i-1,1) = pt(1);
            obj.cp(i-1,2) = pt(2);
        else
            obj.cp(i  ,1) = pt(1);
            obj.cp(i  ,2) = pt(2);
        end
    end
end