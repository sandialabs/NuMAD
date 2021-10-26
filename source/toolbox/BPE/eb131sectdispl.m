function [sectionData] = eb131sectdispl(nodeData,sectionData,flagData,sizedispl)
% sandia equiv beam, calcs the 6 overall displacements at each of the selected sections
% 28 may 2002, mods for interpolating between adjacent fea sections
% version 18, 19 aug 2002, uses a cubic to fit deformed section in z direction
% renamed as eb21sectdispl.m;  includes even powers as well as odd
% 8 may 2003, v 24, weighted solution for section displts using node masses
% 26 june 03, v 25, weighted solution bypassed
% 14 july 2003, wtflag added for control of mass weighting
% 16 july 2003, sizedispl added for control of terms in solution of section displacements
% 17 july 2003, add writing of node displts at section 3
%  20 nov 2003, v27, for use with Legacy  version of Numad
%  14 april 2004, v30, added outer profile
%  22 july 2005, v31, reference to profile removed (solid elements)
% 4 may 2004, corrected temp(3,7) => temp(3,8) etc
% 25 aug 05, test:  added z-displ terms for y^2 x^3 and y^3 x^3
% 2 sept 09, brr: cleaning up the code for use with NuMAD

nbeamnode = nodeData.nbeamnode;
zbeamnode = nodeData.zbeamnode;

sectionnodecount1 = sectionData.sectionnodecount1;
sectionnodecoord1 = sectionData.sectionnodecoord1;
sectionnodedispl1 = sectionData.sectionnodedispl1;
zbeamnode1 = sectionData.zbeamnode1;
sectionmaxchord1 = sectionData.sectionmaxchord1;
sectionchordrot1 = sectionData.sectionchordrot1;

sectionnodecount2 = sectionData.sectionnodecount2;
sectionnodecoord2 = sectionData.sectionnodecoord2;
sectionnodedispl2 = sectionData.sectionnodedispl2;
zbeamnode2 = sectionData.zbeamnode2;
sectionmaxchord2 = sectionData.sectionmaxchord2;
sectionchordrot2 = sectionData.sectionchordrot2;

sectioncofpx1 = sectionData.sectioncofpx1;
sectioncofpx2 = sectionData.sectioncofpx2;
sectioncofpy1 = sectionData.sectioncofpy1;
sectioncofpy2 = sectionData.sectioncofpy2;

sectionLE1 = sectionData.sectionLE1;
sectionLE2 = sectionData.sectionLE2;

sectionnodemass1 = sectionData.sectionnodemass1;
sectionnodemass2 = sectionData.sectionnodemass2;

wtflag = flagData.wtflag;

sectiondispl = zeros(nbeamnode,6,6);
sectionchordrot = zeros(nbeamnode);
sectionmaxchord = zeros(nbeamnode);
sectioncofpx    = zeros(nbeamnode);
sectioncofpy    = zeros(nbeamnode);
sectionLE = zeros(nbeamnode,2);
storesectiondispl = zeros(nbeamnode,1+6*6);
matrixtemp = zeros(sizedispl,6);
matrix1 = zeros(6,6);
matrix2 = zeros(6,6);
matrix  = zeros(6,6);

for isection = 1:nbeamnode % loop on beam sections
    
%     figure(1)
%     plot(sectionnodecoord1(isection,:,2),sectionnodecoord1(isection,:,1),'k*')
%     xlabel('y'); ylabel('x'); axis equal; title(['isection=' num2str(isection)])
%     pause(0.5);
    
    for inode = 1:sectionnodecount1(isection) % loop on fea nodes at near fea section
        temp = zeros(3,sizedispl);
        temp(1,1) = 1;
        temp(2,2) = 1;
        temp(3,3) = 1;
        x = +sectionnodecoord1(isection,inode,1);
        y = +sectionnodecoord1(isection,inode,2);
        temp(1,6) = -y;
        temp(2,6) = +x;
        temp(3,4) = y;
        temp(3,5) = -x;
        if sizedispl > 6 % add  terms as requested  by sizedispl
            temp(3,7) = x*y;
        end
        if sizedispl > 7
            temp(3,8) = x^2;
        end
        if sizedispl > 8
            temp(3,9) = y^2;
        end
        if sizedispl > 9
            temp(3,10) = x^3;
        end
        if sizedispl > 10
            temp(3,11) = y^3;
        end
        if sizedispl > 11
            temp(3,12) = x^2*y;
        end
        if sizedispl > 12
            temp(3,13) = y^2*x;
        end
        if sizedispl > 13
            temp(3,14) = y^3*x^2;
        end
        if sizedispl > 14
            temp(3,15) = x^3*y^2;
        end
        
        if (inode > 1)
            sectionmatrix = [sectionmatrix;temp];  % construct matrix
        else
            sectionmatrix = temp;
        end
        
    end  % end inode loop
    
    % construct weighting matrix
    if wtflag ==1  % add weighting to geom matrix if appropriate
        massvect = [ ];
        for inode = 1:sectionnodecount1(isection)
            vect = ones(3,1)*sectionnodemass1(isection,inode);
            if inode > 1
                massvect = [massvect;vect];
            else
                massvect = vect;
            end
        end   % end of loop on section nodes
        
        weightmatrix = diag(massvect);
        
        %  isection
        %  'weighting matrix, side 1'
        sectionmatrix = weightmatrix*sectionmatrix;
    end %  end wtflag if
    
    for iload = 1:6 % loop on the 6 load cases
        
        for inode = 1:sectionnodecount1(isection) % loop on fea nodes at section
            
            vect = zeros(3,1);  % extract relevant 3 displacements
            for i=1:3
                vect(i,1) = sectionnodedispl1(isection,inode,(iload-1)*3+i);
            end
            
            if (inode > 1)
                nodedisplvect = [nodedisplvect;vect];    % construct vector of displacements
            else
                nodedisplvect = vect;
            end
            
        end % end of inode loop
        
        if iload>1
            nodedisplmatrix = [nodedisplmatrix,nodedisplvect];
        else
            nodedisplmatrix = nodedisplvect;
        end
        
    end  % end of iload loop
    
    if wtflag == 1  % add weighting to vectors of node displts if appropriate
        nodedisplmatrix = weightmatrix* nodedisplmatrix;
    end
    
    matrixtemp = sectionmatrix\nodedisplmatrix;   % divide displ by geom matrix
    matrix1 = matrixtemp(1:6,1:6);                % remove cubic etc term solutions
    
    for inode = 1:sectionnodecount2(isection)      % loop on fea nodes at far side fea section
        temp = zeros(3,sizedispl);
        temp(1,1) = 1;
        temp(2,2) = 1;
        temp(3,3) = 1;
        x = +sectionnodecoord2(isection,inode,1);
        y = +sectionnodecoord2(isection,inode,2);
        temp(1,6) = -y;
        temp(2,6) = +x;
        temp(3,4) = y;
        temp(3,5) = -x;
        if sizedispl > 6   % add  terms as requested  by sizedispl
            temp(3,7) = x*y;
        end
        if sizedispl > 7
            temp(3,8) = x^2;
        end
        if sizedispl > 8
            temp(3,9) = y^2;
        end
        if sizedispl > 9
            temp(3,10) = x^3;
        end
        if sizedispl > 10
            temp(3,11) = y^3;
        end
        if sizedispl > 11
            temp(3,12) = x^2*y;
        end
        if sizedispl > 12
            temp(3,13) = y^2*x;
        end
        if sizedispl > 13
            temp(3,14) = y^3*x^2;
        end
        if sizedispl > 14
            temp(3,15) = x^3*y^2;
        end
        
        if (inode > 1)
            sectionmatrix = [sectionmatrix;temp];  % construct matrix
        else
            sectionmatrix = temp;
        end
        
    end % end inode loop
    
    % construct weighting matrix
    if wtflag ==1   % add weighting to geom matrix if appropriate
        massvect = [];
        for inode = 1:sectionnodecount2(isection)
            vect = ones(3,1)*sectionnodemass2(isection,inode);
            if inode > 1
                massvect = [massvect;vect];
            else
                massvect = vect;
            end
        end % end of loop on section nodes
        
        weightmatrix = diag(massvect);
        
        % isection
        % 'weighting of matrix side 2'
        sectionmatrix = weightmatrix*sectionmatrix;
        
    end  %  end wtflag if
    
    
    for iload = 1:6 % loop on the 6 load cases
        
        for inode = 1:sectionnodecount2(isection)  % loop on fea nodes at section
            
            vect = zeros(3,1); % extract relevant 3 displacements
            for i=1:3
                vect(i,1) = sectionnodedispl2(isection,inode,(iload-1)*3+i);
            end
            
            if (inode > 1)
                nodedisplvect = [nodedisplvect;vect];  % construct vector of displacements
            else
                nodedisplvect = vect;
            end
            
        end % end of inode loop
        
        if iload>1
            nodedisplmatrix = [nodedisplmatrix,nodedisplvect];
        else
            nodedisplmatrix = nodedisplvect;
        end
        
    end                                       % end of iload loop
    
    if wtflag ==1  % add weighting to vectors of node displts if appropriate
        nodedisplmatrix = weightmatrix* nodedisplmatrix;
    end
    
    matrixtemp = sectionmatrix\nodedisplmatrix;   % divide displ by geom matrix
    matrix2 = matrixtemp(1:6,1:6);  % remove cubic term solutions
    
    %interpolate to select displacement at beam node
    dz = zbeamnode2(isection)-zbeamnode1(isection);
    if dz >0
        f = (zbeamnode2(isection)-zbeamnode(isection))/dz;
    else
        f = 0;
    end
    
    matrix = f*matrix1 + (1-f)*matrix2;
    sectiondispl(isection,1:6,1:6) = matrix(:,:);  %put into array
    % interpolate for other  quantities
    sectionchordrot(isection) = f*sectionchordrot1(isection)+(1-f)*sectionchordrot2(isection);
    sectionmaxchord(isection) = f*sectionmaxchord1(isection)+(1-f)*sectionmaxchord2(isection);
    sectioncofpx(isection) = f*sectioncofpx1(isection)+(1-f)*sectioncofpx2(isection);
    sectioncofpy(isection) = f*sectioncofpy1(isection)+(1-f)*sectioncofpy2(isection);
    sectionLE(isection,:) = f*sectionLE1(isection,:) + (1-f)*sectionLE2(isection,:);
    
    for j=1:6    %put into storage in suitable format
        storesectiondispl(isection,(j-1)*6+2:j*6+1) = matrix(1:6,j)';
    end
    
    storesectiondispl(isection,1) = zbeamnode(isection);
    
end % end of isection loop

save sectiondispl.txt storesectiondispl -ascii;

sectionData.sectioncofpx = sectioncofpx;
sectionData.sectioncofpy = sectioncofpy;
sectionData.sectiondispl = sectiondispl;
sectionData.sectionmaxchord = sectionmaxchord;
sectionData.sectionchordrot = sectionchordrot;
sectionData.sectionLE = sectionLE;

end

