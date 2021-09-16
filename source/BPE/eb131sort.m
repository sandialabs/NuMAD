function [sectionData] = eb131sort(feaData,nodeData)
%  sorts the finite element nodes into selected beam node sections
%  modified for random arrangement of fea nodes and double precision in fea values
%
%  28 may 2002, modified for recording both fea sections that are adjacent to specified location
%  latest update: 2 dec 2002, grouping criterion reduced from 0.00001 to 0.0001
%  8 may 2003, v 24, sectionnodemass arrays added for mass weighting of solutions
%  19 aug 2003, detection of near-circle added for max chord and chord orientation
%  20 nov 2003, v27, for use with Legacy Numad
%  latest version: 26 jan 2004, v28
%  latest version: 14 april 2004, v30, added outer profile
%  22 july 2005, dm, v31, solid element version, modifications to remove reference to profile
%  04 aug 2005, dm, change to proximity criterion
%  02 sept 2009, v131, brr, cleaned up for NuMAD usage
%



nfeanode = feaData.nfeanode;
coord    = feaData.coord;
displ    = feaData.displ;
feamass  = feaData.feamass;

nbeamnode = nodeData.nbeamnode;
zbeamnode = nodeData.zbeamnode;


zbeamnode1 = zeros(1,nbeamnode);
zbeamnode2 = zeros(1,nbeamnode);
zsection1 = zeros(1,nbeamnode);
zsection2 = zeros(1,nbeamnode);
sectioncofpx1 = zeros(nbeamnode);
sectioncofpx2 = zeros(nbeamnode);
sectioncofpy1 = zeros(nbeamnode);
sectioncofpy2 = zeros(nbeamnode);
sectionLE1 = zeros(nbeamnode,2);
sectionLE2 = zeros(nbeamnode,2);

% Determine the fea nodes that are adjacent to the spec beam node (on the near and far sides)
for inode = 1:nbeamnode % loop on beam nodes
    
    z = zbeamnode(inode); % z-location of beam node
    ii =1;
    jj = nfeanode;  % number of fea nodes
    dist1 = abs(z); % distance from root to beam node
    dist2 = abs(zbeamnode(nbeamnode)-z); % distance from tip to beam node
    for ifea = 2:nfeanode; % loop on fea nodes
        if z-coord(ifea,3) >= 0; % if fea node is inboard of beam node (inclusive)
            if z-coord(ifea,3) < dist1; 
                ii = ifea;
                dist1 = z-coord(ifea,3);
            end
        end
        if coord(ifea,3)-z >= 0;
            if coord(ifea,3)-z <= dist2
                jj = ifea;
                dist2 = coord(ifea,3)-z;
            end
        end
        
    end %  end of fea node loop
    
    zbeamnode1(inode) = coord(ii,3);
    zbeamnode2(inode) = coord(jj,3);
    
end % end loop on beam nodes

%  scan the fea nodes and place those that fall into selected beamnode
%    locations into the relevant sections

sectionnodecount1 = zeros(nbeamnode,1);  % zero the section counts
sectionnodecount2 = zeros(nbeamnode,1);  % zero the section counts
for ifea = 1:nfeanode % loop on fea nodes
    z = coord(ifea,3);
    for inode=1:nbeamnode % loop on beam nodes
        if abs(z-zbeamnode1(inode))/zbeamnode(nbeamnode)<0.001 % section on near side of beam node, mod 23 sept 05
            sectionnodecount1(inode)=sectionnodecount1(inode)+1;
            sectionnodecoord1(inode,sectionnodecount1(inode),1:3) = coord(ifea,:);
            sectionnodedispl1(inode,sectionnodecount1(inode),1:18) = displ(ifea,:);
            sectionnodemass1(inode,sectionnodecount1(inode)) = feamass(ifea);
            break;
        end
    end  % end loop on beam nodes
    
    for inode=1:nbeamnode  % loop on beam nodes
        if abs(z-zbeamnode1(inode))/zbeamnode(nbeamnode)<0.001 % section on near side of beam node, mod 23 sept 05
            sectionnodecount2(inode)=sectionnodecount2(inode)+1;
            sectionnodecoord2(inode,sectionnodecount2(inode),1:3) = coord(ifea,:);
            sectionnodedispl2(inode,sectionnodecount2(inode),1:18) = displ(ifea,:);
            sectionnodemass2(inode,sectionnodecount2(inode)) = feamass(ifea);
            break;
        end
    end % end loop on beam nodes
end % end loop on fea nodes

% check each section for chord length and orientation
for inode=1:nbeamnode % loop on beam nodes
    
    maxlength= zeros(1,sectionnodecount1(inode));
    for jnode=1:sectionnodecount1(inode) %loop on section nodes
        x1 = sectionnodecoord1(inode,jnode,1);
        y1 = sectionnodecoord1(inode,jnode,2);
        for knode=1:sectionnodecount1(inode) %loop on section nodes, first of section pair
            x2 = sectionnodecoord1(inode,knode,1);
            y2 = sectionnodecoord1(inode,knode,2);
            length = sqrt((x2-x1)^2+(y2-y1)^2);
            if length > maxlength(jnode)
                maxlength(jnode) = length;
            end
            
        end  % end of knode loop
        
    end % end of jnode loop
    
    [nr,nc] = size(maxlength);
    [maxmaxlength,maxi] = max(maxlength);
    minmaxlength = min(maxlength);
    maxlength = sort(maxlength); % sorts in ascending order
    
    if numel(maxlength) < 12
        error(sprintf('BPE ERROR: Segment edge #%i doesn''t contain enough FE nodes for BPE analysis.  Remesh model with smaller elements or redefine BPE segment edges to avoid this section.',inode));
    end
    
    if maxmaxlength/maxlength(nc-10) < 1.01 % check for near-circle
        disp('near circle')
        xx1=0.0 ;
        xx2=0.0;
        yy1 = min(sectionnodecoord1(inode,:,2));
        yy2 = max(sectionnodecoord1(inode,:,2));
    else
        x1 = sectionnodecoord1(inode,maxi,1);
        y1 = sectionnodecoord1(inode,maxi,2);
        length = 0.0;
        for knode=1:sectionnodecount1(inode) %loop on section nodes, first of section pair
            x2 = sectionnodecoord1(inode,knode,1);
            y2 = sectionnodecoord1(inode,knode,2);
            len = sqrt((x2-x1)^2+(y2-y1)^2);
            if len > length
                length = len;
                kk = knode;
                xx1 = x1;
                xx2 = x2;
                yy1 = y1;
                yy2 = y2;
            end
            
        end % end of knode loop
        
    end % end of if maxlength/....
    
    if yy1 > yy2; %   make xx1, yy1 = leading edge
        temp = xx1;
        xx1 = xx2;
        xx2 = temp;
        temp = yy1;
        yy1 = yy2;
        yy2 = temp;
    end   % end of leading edge check
    
    if abs(yy2-yy1) == 0
        sectionchordrot1(inode) = 0;
    else
        sectionchordrot1(inode) = -atan((xx1-xx2)/(yy1-yy2));   % rotation (right hand rule)
    end
    sectionmaxchord1(inode) = maxmaxlength;
    sectioncofpx1(inode) = xx1 - (xx1-xx2)/4;
    sectioncofpy1(inode) = yy1 - (yy1-yy2)/4;
    
    sectionLE1(inode,1) = xx1;  % store leading edge coords
    sectionLE1(inode,2) = yy1;
    
    % repeat for second of section pairs
    maxlength= zeros(1,sectionnodecount2(inode));
    for jnode=1:sectionnodecount2(inode)          %loop on section nodes
        x1 = sectionnodecoord2(inode,jnode,1);
        y1 = sectionnodecoord2(inode,jnode,2);
        for knode=1:sectionnodecount2(inode)       %loop on section nodes, first of section pair
            x2 = sectionnodecoord2(inode,knode,1);
            y2 = sectionnodecoord2(inode,knode,2);
            length = sqrt((x2-x1)^2+(y2-y1)^2);
            if length > maxlength(jnode)
                maxlength(jnode) = length;
            end
            
        end % end of knode loop
        
    end  % end of jnode loop
    
    [nr,nc] = size(maxlength);
    [maxmaxlength,maxi] = max(maxlength);
    minmaxlength = min(maxlength);
    maxlength = sort(maxlength);   % sorts in ascending order
    if maxmaxlength/maxlength(nc-10) < 1.01   % check for near-circle
        %    'near circle'
        xx1=0.0 ;
        xx2=0.0;
        yy1 = min(sectionnodecoord2(inode,:,2));
        yy2 = max(sectionnodecoord2(inode,:,2));
    else
        x1 = sectionnodecoord2(inode,maxi,1);
        y1 = sectionnodecoord2(inode,maxi,2);
        length = 0.0;
        for knode=1:sectionnodecount2(inode)  %loop on section nodes, first of section pair
            x2 = sectionnodecoord2(inode,knode,1);
            y2 = sectionnodecoord2(inode,knode,2);
            len = sqrt((x2-x1)^2+(y2-y1)^2);
            if len > length
                length = len;
                kk = knode;
                xx1 = x1;
                xx2 = x2;
                yy1 = y1;
                yy2 = y2;
            end
            
        end  % end of knode loop
        
    end   % end of if maxlength/....
    
    
    if yy1 > yy2;   %   make xx1, yy1 = leading edge
        temp = xx1;
        xx1 = xx2;
        xx2 = temp;
        temp = yy1;
        yy1 = yy2;
        yy2 = temp;
    end
    
    if abs(yy2-yy1) == 0
        sectionchordrot2(inode) = 0;
    else
        sectionchordrot2(inode) = -atan((xx1-xx2)/(yy1-yy2));   % rotTwist (right hand rule)
    end
    sectionmaxchord2(inode) = maxmaxlength;
    sectioncofpx2(inode) = xx1 - (xx1-xx2)/4;
    sectioncofpy2(inode) = yy1 - (yy1-yy2)/4;
    
    sectionLE2(inode,1) = xx1;  % store leading edge coords
    sectionLE2(inode,2) = yy1;
end   % end of inode loop

sectionData.sectionnodecount1 = sectionnodecount1;
sectionData.sectionnodecoord1 = sectionnodecoord1;
sectionData.sectionnodedispl1 = sectionnodedispl1;
sectionData.zbeamnode1 = zbeamnode1;
sectionData.sectionmaxchord1 = sectionmaxchord1;
sectionData.sectionchordrot1 = sectionchordrot1;

sectionData.sectionnodecount2 = sectionnodecount2;
sectionData.sectionnodecoord2 = sectionnodecoord2;
sectionData.sectionnodedispl2 = sectionnodedispl2;
sectionData.zbeamnode2 = zbeamnode2;
sectionData.sectionmaxchord2 = sectionmaxchord2;
sectionData.sectionchordrot2 = sectionchordrot2;

sectionData.sectioncofpx1 = sectioncofpx1;
sectionData.sectioncofpx2 = sectioncofpx2;
sectionData.sectioncofpy1 = sectioncofpy1;
sectionData.sectioncofpy2 = sectioncofpy2;

sectionData.sectionLE1 = sectionLE1;
sectionData.sectionLE2 = sectionLE2;

sectionData.sectionnodemass1 = sectionnodemass1;
sectionData.sectionnodemass2 = sectionnodemass2;


end
