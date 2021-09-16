function [storeData] = eb131singular(feaData,nodeData,sectionData,flagData,singularData,storeData,sizedispl)
%
%  lengthens each singular element until non-singular
%  latest version: 26 jan 2004, v28, dm,
%  latest version: 14 april 2004, v30, dm, added outer section
%  latest version: 27 april 2004, v30, dm, new K returned to storestiffsym
% 4 may 2004, corrected temp(3,7) => temp(3,8) etc
%  22 july 2005, v32, reference to profile removed (solid elements)

nbeamnode = nodeData.nbeamnode;
zbeamnode = nodeData.zbeamnode;
wtflag = flagData.wtflag;

sectiondispl = sectionData.sectiondispl;

nsingular = singularData.nsingular;
singular = singularData.singular;
fid6 = singularData.fid6;

nfeanode = feaData.nfeanode;
coord    = feaData.coord;
displ    = feaData.displ;
feamass  = feaData.feamass;

maxsingularcount = 6;       % set max no of length increases

%                        % loop on singular elements
%
for isingular = 1:nsingular
    disp(sprintf('Attempting to deal with singular element #%i of %i...',isingular,nsingular))

    ielement = singular(isingular);
    singularflag = 1;                    %set flag on
    singularcount = 0;
    
    znode1 = zbeamnode(ielement);
    znode2 = zbeamnode(ielement+1);
    length = znode2-znode1;
    origlength = length;                  % set original length
    while singularflag ==1             % check singularity
        singularcount = singularcount + 1;
        if singularcount > maxsingularcount,
            break;
        end
        znode1 = znode1-length*0.25;
        znode2 = znode2 + length*0.25;
        if znode1<0               % check for start of model
            znode1 = 0.0;
        elseif znode2 > zbeamnode(nbeamnode)   % check for end of model
            znode2 = zbeamnode(nbeamnode);
        end
        length = znode2 - znode1;
        
        %  determine the fea nodes that are adjacent to the spec beam node (on the near and far sides)
        %
        for inode = 1:2      % loop on beam nodes
            %
            if inode==1,
                z = znode1;
            elseif inode==2,
                z = znode2;
            end
            ii =1;
            jj = nfeanode;
            dist1 = abs(z);
            %       inode
            dist2 = abs(zbeamnode(nbeamnode)-z);
            for ifea = 2:nfeanode;               % loop on fea nodes
                if z-coord(ifea,3) >=0;
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
                %
            end                           %  end of fea node loop
            %
            %       inode
            z1(inode) = coord(ii,3);       % store z locations of fea
            z2(inode) = coord(jj,3);
            %
            %
        end                             % end loop on beam nodes
        %
        %      scan the fea nodes and place those that fall into selected beamnode
        %        locations into the relevant sections
        %
        singsectionnodecount1 = zeros(2,1);  % zero the section counts
        singsectionnodecount2 = zeros(2,1);  % zero the section counts
        for ifea = 1:nfeanode                                  % loop on fea nodes
            z = coord(ifea,3);
            for inode=1:2                % loop on beam nodes (first of sections)
                if abs(z-z1(inode))/zbeamnode(nbeamnode)<0.0001              % section on near side of beam node, mod 26 jan 2004
                    singsectionnodecount1(inode)=singsectionnodecount1(inode)+1;
                    singsectionnodecoord1(inode,singsectionnodecount1(inode),1:3) = coord(ifea,:);
                    singsectionnodedispl1(inode,singsectionnodecount1(inode),1:18) = displ(ifea,:);
                    singsectionnodemass1(inode,singsectionnodecount1(inode)) = feamass(ifea);
                    break;
                end
            end                                                  % end loop on beam nodes
            
            for inode=1:2              % loop on beam nodes (second of sections)
                %         if abs(z-z2(inode))<0.00001                                     % section on far side of beam node
                if abs(z-z2(inode))/zbeamnode(nbeamnode)<0.0001              % section on far side of beam node,mod 26 jan 2004
                    %        ifea
                    singsectionnodecount2(inode)=singsectionnodecount2(inode)+1;
                    singsectionnodecoord2(inode,singsectionnodecount2(inode),1:3) = coord(ifea,:);
                    singsectionnodedispl2(inode,singsectionnodecount2(inode),1:18) = displ(ifea,:);
                    singsectionnodemass2(inode,singsectionnodecount2(inode)) = feamass(ifea);
                    break
                end
            end                                                  % end loop on beam nodes
            %
            %
        end                             % end loop on fea nodes
        
        
        %                                        % calc new section displacements for this modified beam
        
        sectiondispl = zeros(nbeamnode,6,6);
        storesectiondispl = zeros(nbeamnode,1+6*6);
        matrixtemp = zeros(sizedispl,6);
        matrix1 = zeros(6,6);
        matrix2 = zeros(6,6);
        matrix  = zeros(6,6);
        
        %
        for isection = 1:2                  % loop on start, end of beam
            %
            for inode = 1:singsectionnodecount1(isection)      % loop on fea nodes at near fea section
                temp = zeros(3,sizedispl);
                temp(1,1) = 1;
                temp(2,2) = 1;
                temp(3,3) = 1;
                x = +singsectionnodecoord1(isection,inode,1);
                y = +singsectionnodecoord1(isection,inode,2);
                temp(1,6) = -y;
                temp(2,6) = +x;
                temp(3,4) = y;
                temp(3,5) = -x;
                if sizedispl > 6                 % add  terms as requested  by sizedispl
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
                
                
                if (inode > 1)
                    sectionmatrix = [sectionmatrix;temp];  % construct matrix
                else
                    sectionmatrix = temp;
                end
                %
            end                                               % end inode loop
            %
            %                 construct weighting matrix
            if wtflag ==1                     % add weighting to geom matrix if appropriate
                massvect = [ ];
                for inode = 1:singsectionnodecount1(isection)
                    vect = ones(3,1)*singsectionnodemass1(isection,inode);
                    if inode > 1
                        massvect = [massvect;vect];
                    else
                        massvect = vect;
                    end
                end                              % end of loop on section nodes
                
                weightmatrix = diag(massvect);
                sectionmatrix = weightmatrix*sectionmatrix;
            end                               %  end wtflag if
            
            %
            for iload = 1:6                        % loop on the 6 load cases
                %
                for inode = 1:singsectionnodecount1(isection)    % loop on fea nodes at section
                    %
                    vect = zeros(3,1);                             % extract relevant 3 displacements
                    for i=1:3
                        vect(i,1) = singsectionnodedispl1(isection,inode,(iload-1)*3+i);
                    end
                    %
                    if (inode > 1)
                        nodedisplvect = [nodedisplvect;vect];    % construct vector of displacements
                    else
                        nodedisplvect = vect;
                    end
                    %
                end                                  % end of inode loop
                
                if iload>1
                    nodedisplmatrix = [nodedisplmatrix,nodedisplvect];
                else
                    nodedisplmatrix = nodedisplvect;
                end
                %
            end                                             % end of iload loop
            
            
            %
            if wtflag ==1                     % add weighting to vectors of node displts if appropriate
                nodedisplmatrix = weightmatrix* nodedisplmatrix;
            end
            
            matrixtemp = sectionmatrix\nodedisplmatrix;   % divide displ by geom matrix
            matrix1 = matrixtemp(1:6,1:6);                % remove cubic etc term solutions
            %
            %                                       % repeat for far side section
            %
            for inode = 1:singsectionnodecount2(isection)      % loop on fea nodes at far side fea section
                temp = zeros(3,sizedispl);
                temp(1,1) = 1;
                temp(2,2) = 1;
                temp(3,3) = 1;
                x = +singsectionnodecoord2(isection,inode,1);
                y = +singsectionnodecoord2(isection,inode,2);
                temp(1,6) = -y;
                temp(2,6) = +x;
                temp(3,4) = y;
                temp(3,5) = -x;
                if sizedispl > 6                 % add  terms as requested  by sizedispl
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
                
                if (inode > 1)
                    sectionmatrix = [sectionmatrix;temp];  % construct matrix
                else
                    sectionmatrix = temp;
                end
                %
            end                                             % end inode loop
            
            %                   construct weighting matrix
            if wtflag ==1                     % add weighting to geom matrix if appropriate
                massvect = [];
                for inode = 1:singsectionnodecount2(isection)
                    vect = ones(3,1)*singsectionnodemass2(isection,inode);
                    if inode > 1
                        massvect = [massvect;vect];
                    else
                        massvect = vect;
                    end
                end                              % end of loop on inode / section nodes
                
                weightmatrix = diag(massvect);
                sectionmatrix = weightmatrix*sectionmatrix;
                
            end                              %  end wtflag if
            
            %
            for iload = 1:6                        % loop on the 6 load cases
                %
                for inode = 1:singsectionnodecount2(isection)    % loop on fea nodes at section
                    %
                    vect = zeros(3,1);                             % extract relevant 3 displacements
                    for i=1:3
                        vect(i,1) = singsectionnodedispl2(isection,inode,(iload-1)*3+i);
                    end
                    %
                    if (inode > 1)
                        nodedisplvect = [nodedisplvect;vect];    % construct vector of displacements
                    else
                        nodedisplvect = vect;
                    end
                    %
                end                                  % end of inode loop
                
                if iload>1
                    nodedisplmatrix = [nodedisplmatrix,nodedisplvect];
                else
                    nodedisplmatrix = nodedisplvect;
                end
                %
            end                                             % end of iload loop
            
            %
            if wtflag ==1                     % add weighting to vectors of node displts if appropriate
                nodedisplmatrix = weightmatrix* nodedisplmatrix;
            end
            
            matrixtemp = sectionmatrix\nodedisplmatrix;   % divide displ by geom matrix
            matrix2 = matrixtemp(1:6,1:6);        % remove cubic term solutions
            %
            %                                            %interpolate to select displacement at beam node
            dz = sectionData.zbeamnode2(isection)-sectionData.zbeamnode1(isection);
            if dz >0
                f = (sectionData.zbeamnode2(isection)-sectionData.zbeamnode(isection))/dz;
            else
                f = 0;
            end
            %
            matrix = f*matrix1 + (1-f)*matrix2;
            sectiondispl(isection,1:6,1:6) = matrix(:,:);  %put into array
            
        end                                    % end of isection=1:2 loop
        
        %                                              % start of stiffness matrix formulation
        %                                            % create matrix of relative displacements
        for iload = 1:6                         % loop on load cases
            %
            relativedispl(1,iload) = sectiondispl(2,1,iload)...
                -sectiondispl(1,1,iload)-length*sectiondispl(1,5,iload);
            relativedispl(2,iload) = sectiondispl(2,2,iload)...
                -sectiondispl(1,2,iload)+length*sectiondispl(1,4,iload);
            relativedispl(3,iload) = sectiondispl(2,3,iload)...
                -sectiondispl(1,3,iload);
            relativedispl(4,iload) = sectiondispl(2,4,iload)...
                -sectiondispl(1,4,iload);
            relativedispl(5,iload) = sectiondispl(2,5,iload)...
                -sectiondispl(1,5,iload);
            relativedispl(6,iload) = sectiondispl(2,6,iload)...
                -sectiondispl(1,6,iload);
            %
        end                              % end of loop on load cases
        %
        %                               % create matrix of loads at end 2 of element
        %
        force = zeros(6,6);
        for iload = 1:6           % loop on loads
            %
            force(iload,iload) = 1;
            if iload == 1
                force(5,iload) = zbeamnode(nbeamnode) - znode2;
            end
            if iload == 2
                force(4,iload) = -zbeamnode(nbeamnode) + znode2;
            end
            %
        end                              % end of loop on loads
        %
        %          % create stiffness matrix by dividing forces by displacements
        %
        matrix = force/relativedispl;
        matrix = (matrix + matrix')/2;          % make matrix symmetric
        %                                            % check matrix for singularities
        try
            lambda = eig(matrix);
        catch
            error(['ERROR....Choose new BPE nodes.  Problem in BPE element #' num2str(ielement)])
            break;
        end
        
        mineig = min(lambda);
        if mineig > 0
            singularflag = 0;
            %                                            % convert to section stiffness
            % %reate common E matrix
            E = zeros(6,6);
            E(1,5) = 1;
            E(2,4) = -1;
            
            vector = [1 1 1 1 1 1]*length;
            %   % create H matrix
            H = diag(vector);
            H(4,2) = -(length^2)/2;
            H(5,1) = -H(4,2);
            %  % create L matrix
            L = diag(vector)*length/2;
            L(4,2) = -(length^3)/3;
            L(5,1) = -L(4,2);
            %  % inputs for Lyapunov eqtn: -C = AX + XB
            B = H*inv(L);
            A = E;
            C = matrix;
            C = inv(C)*inv(L);
            X = lyap(A,B,C);                           %  solve Lyapunov equation for X =-inv(stiffness matrix)
            matrix = -inv(X);
            %
            %                                                   %now generate element stiffness for original length
            vector = [1 1 1 1 1 1]*origlength;      % create H matrix
            H = diag(vector);
            H(4,2) = -(origlength^2)/2;
            H(5,1) = - H(4,2);
            Q = diag(vector)*origlength/2;          % create Q matrix
            Q(4,2) = -(origlength^3)/3;
            Q(5,1) = -Q(4,2);
            matrix = inv(inv(matrix)*H + E*inv(matrix)*Q);
            
            i1 =(ielement-1)*8;
           storeData.storestiffsym(i1+1:i1+6,1:6) = matrix;          % put new K back in storestiffsym
            
        end                   % end of if check for singularity
        
        fprintf(fid6,'element # %2i cycle count # %2i length = %9.3f\n',ielement,singularcount,length);
        
        if singularflag==1                                 % print status of singularity
            fprintf(fid6,'element is still singular\n');
            disp(sprintf('End of iteration #%i...Element still singular at singular element #%i of %i.',singularcount,isingular,nsingular))
        else
            fprintf(fid6,'element is no longer singular\n');
            disp(sprintf('Element is no longer singular at singular element #%i of %i.',isingular,nsingular))
        end
        
    end                     % end of while loop
    
end                      % end of loop on singular elements

% storeData.storestiffsym = storestiffsym;

end


