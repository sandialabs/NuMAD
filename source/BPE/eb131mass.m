function [massData] = eb131mass(feaData,nodeData,sectionData,flagData,kData)
%  creates all required mass and inertia information
%  latest update  11 July 2003, massmts about centroid, cross products = 0
%  20 nov 2003, v27, for use with Legacy Numad
%  13 april 2004, v30, profiles
%  22 july 2005, v31, reference to profile changed back to section

nfeanode = feaData.nfeanode;
feamass = feaData.feamass;
coord   = feaData.coord;

nbeamnode = nodeData.nbeamnode;
zbeamnode = nodeData.zbeamnode;

sectioncofpx = sectionData.sectioncofpx;
sectioncofpy = sectionData.sectioncofpy;
sectionLE = sectionData.sectionLE;
sectionmaxchord = sectionData.sectionmaxchord;
sectionchordrot = sectionData.sectionchordrot;

feaflag = flagData.feaflag;

givenk = kData.givenk;
npairk = kData.npairk;
zpairk = kData.zpairk;

nodemass = zeros(nbeamnode,1);
massmtx = zeros(nbeamnode,1);
massmty = zeros(nbeamnode,1);
massmtz = zeros(nbeamnode,1);
massmtx2 = zeros(nbeamnode,1);
massmty2 = zeros(nbeamnode,1);
massmtz2 = zeros(nbeamnode,1);
massmtxy = zeros(nbeamnode,1);
massmtxz = zeros(nbeamnode,1);
massmtyz = zeros(nbeamnode,1);
Ixx = zeros(nbeamnode,1);
Iyy = zeros(nbeamnode,1);
Ixy = zeros(nbeamnode,1);
cgx = zeros(nbeamnode,1);
cgy = zeros(nbeamnode,1);
cgz = zeros(nbeamnode,1);
%  the fea model is scanned and each inertia property is added to appropriate element
%
if feaflag==1.0                % check if fea has been read
    for ifea = 1:nfeanode                     %  loop on all fea nodes
        
        zfea = coord(ifea,3);                          % spanwise coordinate
        if zfea > 0                         % neglect all nodes at root (or negative)
            inode = nbeamnode-1;               % start at last but one node
            while zfea <= zbeamnode(inode)     % find beam element
                inode = inode-1;
            end                                             % end of while loop
            %
            length = zbeamnode(inode+1) - zbeamnode(inode);
            nodemass(inode) = nodemass(inode) + feamass(ifea);               % total mass of element
            % inertia tensor
            massmtx(inode) = massmtx(inode) + coord(ifea,1)*feamass(ifea);
            massmty(inode) = massmty(inode) + coord(ifea,2)*feamass(ifea);
            massmtz(inode) = massmtz(inode) + (coord(ifea,3)-zbeamnode(inode))*feamass(ifea);
            massmtz2(inode) = massmtz2(inode) + (coord(ifea,1)^2+coord(ifea,2)^2)*feamass(ifea);
            massmty2(inode) = massmty2(inode) + ((coord(ifea,3)-zbeamnode(inode))^2+coord(ifea,1)^2)*feamass(ifea);
            massmtx2(inode) = massmtx2(inode) + ((coord(ifea,3)-zbeamnode(inode))^2+coord(ifea,2)^2)*feamass(ifea);
            
            massmtxy(inode) = massmtxy(inode) + coord(ifea,1)*coord(ifea,2)*feamass(ifea);
            massmtxz(inode) = massmtxz(inode) + coord(ifea,1)*(coord(ifea,3)-zbeamnode(inode))*feamass(ifea);
            massmtyz(inode) = massmtyz(inode) + coord(ifea,2)*(coord(ifea,3)-zbeamnode(inode))*feamass(ifea);
            Ixx(inode) = Ixx(inode) + (coord(ifea,2)^2*feamass(ifea))/length;
            Iyy(inode) = Iyy(inode) + (coord(ifea,1)^2*feamass(ifea))/length;
            Ixy(inode) = Ixy(inode) + (coord(ifea,1)*coord(ifea,2)*feamass(ifea))/length;
            %
            
        end;                               % end of if check for z > 0
        
    end;                               % end of loop on fea nodes
    %                                 % calc center of mass of each element
    for inode = 1:nbeamnode-1           % loop on elements
        
        %
        cgx(inode) = massmtx(inode)/nodemass(inode);
        cgy(inode) = massmty(inode)/nodemass(inode);
        cgz(inode) = massmtz(inode)/nodemass(inode);
        
        %  adjust mass mts to be about c of g
        
        massmtx2(inode) = massmtx2(inode) - nodemass(inode)*cgz(inode)^2;
        massmty2(inode) = massmty2(inode) - nodemass(inode)*cgz(inode)^2;
        
        if massmtx2(inode)<0            % check for negative inertias
            disp('negative Ixx')
        elseif massmty2(inode)<0
            disp('negative Iyy')
        end
        
        massmtxy(inode) =0;             %  put cross products = 0
        massmtxz(inode) =0;
        massmtyz(inode) =0;
        
    end;                               % end of loops on elements
    
else         % if no fea model exists
    
    %             now check if elements are within specified regions
    for ielement = 1:nbeamnode-1         % loop on elements
        %
        flaggivenk1 = 0;                  % check for use of given k
        for ipairk = 1:npairk    % check start of beam
            if zbeamnode(ielement) >= zpairk(ipairk,1)
                if zbeamnode(ielement) <= zpairk(ipairk,2)
                    flaggivenk1 = 1;
                    ipair1 = ipairk;
                    break;
                end
            end                         % end of check for given k
        end                         % end of ipairk loop
        %
        flaggivenk2 = 0;
        for ipairk = 1:npairk    % check end of beam
            if zbeamnode(ielement+1) >= zpairk(ipairk,1)
                if zbeamnode(ielement+1) <= zpairk(ipairk,2)
                    flaggivenk2 = 1;
                    ipair2 = ipairk;
                    break;
                end
            end                         % end of check for given k
        end                         % end of ipairk loop
        %
        if flaggivenk1 == 1 && flaggivenk2==1
            z = zbeamnode(ielement);        % extract section inertias per len for start of beam element
            f = (z - zpairk(ipair1,1))/(zpairk(ipair1,2)-zpairk(ipair1,1));   % fraction for proportioning
            
            nodemass(ielement) = (1-f)*givenk((ipair1-1)*16+1,2) + f*givenk((ipair1-1)*16+9,2);
            sectioncofpy(ielement) =  (1-f)*givenk((ipair1-1)*16+1,3) + f*givenk((ipair1-1)*16+9,3) ;
            
            Ixx(ielement) =  (1-f)*givenk((ipair1-1)*16+1,4) + f*givenk((ipair1-1)*16+9,4) ;
            Iyy(ielement) =  (1-f)*givenk((ipair1-1)*16+1,5) + f*givenk((ipair1-1)*16+9,5) ;
            Izz(ielement) =  (1-f)*givenk((ipair1-1)*16+1,6) + f*givenk((ipair1-1)*16+9,6) ;
            
            cgx(ielement) = (1-f)*givenk((ipair1-1)*16+2,1) + f*givenk((ipair1-1)*16+10,1);
            cgy(ielement) = (1-f)*givenk((ipair1-1)*16+2,2) + f*givenk((ipair1-1)*16+10,2);
            sectionmaxchord(ielement) = (1-f)*givenk((ipair1-1)*16+2,3) + f*givenk((ipair1-1)*16+10,3);
            sectionchordrot(ielement) = (1-f)*givenk((ipair1-1)*16+2,4) + f*givenk((ipair1-1)*16+10,4);
            sectionLE(ielement,1) = (1-f)*givenk((ipair1-1)*16+2,5) + f*givenk((ipair1-1)*16+10,5);
            sectionLE(ielement,2) = (1-f)*givenk((ipair1-1)*16+2,6) + f*givenk((ipair1-1)*16+10,6);
            
            %
            z = zbeamnode(ielement+1);                                      % extract section inertias for end of beam element
            f = (z - zpairk(ipair2,1))/(zpairk(ipair2,2)-zpairk(ipair2,1));   % fraction for proportioning
            
            %                                                                       % add to the first and  average
            nodemass(ielement) = (nodemass(ielement)+(1-f)*givenk((ipair2-1)*16+1,2) + f*givenk((ipair2-1)*16+9,2))/2;
            sectioncofpy(ielement) = (sectioncofpy(ielement)+(1-f)*givenk((ipair2-1)*16+1,3) + f*givenk((ipair2-1)*16+9,3))/2;
            Ixx(ielement) = (Ixx(ielement)+(1-f)*givenk((ipair2-1)*16+1,4) + f*givenk((ipair2-1)*16+9,4))/2;
            Iyy(ielement) = (Iyy(ielement)+(1-f)*givenk((ipair2-1)*16+1,5) + f*givenk((ipair2-1)*16+9,5))/2;
            Izz(ielement) = (Izz(ielement)+(1-f)*givenk((ipair2-1)*16+1,6) + f*givenk((ipair2-1)*16+9,6))/2;
            cgx(ielement) = (cgx(ielement)+(1-f)*givenk((ipair2-1)*16+2,1) + f*givenk((ipair2-1)*16+10,1))/2;
            cgy(ielement) = (cgy(ielement)+(1-f)*givenk((ipair2-1)*16+2,2) + f*givenk((ipair2-1)*16+10,2))/2;
            sectionmaxchord(ielement) = (sectionmaxchord(ielement)+(1-f)*givenk((ipair2-1)*16+2,3) + f*givenk((ipair2-1)*16+10,3))/2;
            sectionchordrot(ielement) = (sectionchordrot(ielement)+(1-f)*givenk((ipair2-1)*16+2,4) + f*givenk((ipair2-1)*16+10,4))/2;
            sectionLE(ielement,1) = (sectionLE(ielement,1)+(1-f)*givenk((ipair2-1)*16+2,5) + f*givenk((ipair2-1)*16+10,5))/2;
            sectionLE(ielement,2) = (sectionLE(ielement,2)+(1-f)*givenk((ipair2-1)*16+2,6) + f*givenk((ipair2-1)*16+10,6))/2;
            
            %    now adjust values from per unit length values to values for element, assuming uniform properties, about element CofG
            length = zbeamnode(ielement+1) - zbeamnode(ielement);
            
            cgz(ielement) = length/2;
            nodemass(ielement) = nodemass(ielement)*length;
            massmtx2(ielement) = Ixx(ielement)*length + nodemass(ielement)*length^2/12;       % mod 11 july 2003
            massmty2(ielement) = Iyy(ielement)*length + nodemass(ielement)*length^2/12;       % mod 11 july 2003
            massmtz2(ielement) = Izz(ielement)*length ;
            massmtxy(ielement) = 0 ;
            massmtxz(ielement) = 0 ;
            massmtyz(ielement) = 0 ;
            sectioncofpx(ielement) = 0;
            
        end                                  % end of if for flags for givenk
        
    end                                  %  end of loop on elements
    
end                                % end of if.. else for fea flag


massData.nodemass = nodemass;
massData.cgx = cgx;
massData.cgy = cgy;
massData.cgz = cgz;
massData.massmtx2 = massmtx2;
massData.massmty2 = massmty2;
massData.massmtz2 = massmtz2;
massData.massmtxy = massmtxy;
massData.massmtxz = massmtxz;
massData.massmtyz = massmtyz;
massData.Ixx = Ixx;
massData.Iyy = Iyy;
massData.Ixy = Ixy;

end
