function midNodei=findCenterNode(nodes,direction)        
    xlist=nodes(:,direction+1); %list of nodal x values
    [xmax,~]=max(xlist);
    [xmin,~]=min(xlist);

    xmid=mean([xmin,xmax]); %x coordinate of sparcap centerline

    nnode=numel(nodes(:,1));
    minDiff=100; %initialize
    for iNode=1:nnode
        x=nodes(iNode,direction+1);    %x position of node
        diff = abs(x-xmid); %distance between the node and sparcap centerline
        if diff < minDiff
          midNodei=iNode;
          minDiff =diff;
        end

    end
end