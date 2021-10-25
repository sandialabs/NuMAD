function [storeData,singularData] = eb131stiff(nodeData,sectionData,flagData,kData)
% sandia equiv beam study
% calculates the 6x6 stiffnesses of each beam element (end 2 wrt end 1)
% version 20, 10 oct 02: over-ride by specified k matrices where given
% latest: 27 feb 2003, check both ends of beam for givenk pairs
% latest: 14july 2003, expflag check added for exponential interpolation
% 20 nov 2003, v27, for use with legacy Numad
% 16 jan 2004, v28, for use with singular routine
% 22 jan 2004, changed "displ" to "relativedispl" to avoid overwriting initial global displ matrix
% 26 jan 2004, v28, add some print statements
%  14 april 2004, v30, added outer profile
% 24 sept  2004, v30, switched k1 and k2 in definitions of kk1 & kk2
%  22 july 2005, v31, reference to profile removed (solid elements)
% 02 sept 2009, brr, cleaned up for NuMAD integration
%

givenk = kData.givenk;
npairk = kData.npairk;
zpairk = kData.zpairk;

singular = [];

feaflag = flagData.feaflag;
expflag = flagData.expflag;

nbeamnode = nodeData.nbeamnode;
zbeamnode = nodeData.zbeamnode;

sectiondispl = sectionData.sectiondispl;

storestiff = zeros((nbeamnode)*8,6);
storestiffsym = zeros((nbeamnode)*8,6);
storeflexsym = zeros((nbeamnode)*8,6);
eigstore = zeros(6,nbeamnode-1);
matrix = zeros(6,6);
flex = zeros(6,6);
E = zeros(6,6);  % common E matrix
E(1,5) = 1;
E(2,4) = -1;
nsingular = 0;
fid6 = fopen('singularout.txt','w');

for ielement = 1:nbeamnode-1       % loop on elements
    lambda = zeros(6,1);
    relativedispl = zeros(6,6);
    length = zbeamnode(ielement+1)-zbeamnode(ielement);
    
    flaggivenk1 = 0;  % check for use of given k
    for ipairk = 1:npairk    % check start of beam
        if zbeamnode(ielement) >= zpairk(ipairk,1)
            if zbeamnode(ielement) <= zpairk(ipairk,2)
                flaggivenk1 = 1;
                ipair1 = ipairk;
                break;
            end
        end   % end of check for given k
    end   % end of ipairk loop
    
    flaggivenk2 = 0;
    for ipairk = 1:npairk    % check end of beam
        if zbeamnode(ielement+1) >= zpairk(ipairk,1)
            if zbeamnode(ielement+1) <= zpairk(ipairk,2)
                flaggivenk2 = 1;
                ipair2 = ipairk;
                break;
            end
        end   % end of check for given k
    end    % end of ipairk loop
    
    if flaggivenk1 == 1 && flaggivenk2==1
        z = zbeamnode(ielement);
        
        z1 = zpairk(ipair1,1);          % calc k for near end of beam
        z2 = zpairk(ipair1,2);
        f = (z-z1)/(z2-z1);
        k1 = givenk((ipair1-1)*16+3:(ipair1-1)*16+8,:);
        k2 = givenk((ipair1-1)*16+11:(ipair1-1)*16+16,:);
        
        if expflag ==1                        % use exponential interpolation
            for i=1:6
                for j=1:6
                    kk1(i,j) = k1(i,j)*(k2(i,j)/k1(i,j))^f ;
                end
            end
        else                                   % use linear interpolation
            kk1 = f*k2 + (1-f)*k1;             % 24 sept 04, k1, k2 switched
        end
        
        z = zbeamnode(ielement+1);
        z1 = zpairk(ipair2,1);          % calc k for far end of beam
        z2 = zpairk(ipair2,2);
        f = (z-z1)/(z2-z1);
        k1 = givenk((ipair2-1)*16+3:(ipair2-1)*16+8,:);
        k2 = givenk((ipair2-1)*16+11:(ipair2-1)*16+16,:);
        
        if expflag ==1                        % use exponential interpolation
            for i=1:6
                for j=1:6
                    kk2(i,j) = k1(i,j)*(k2(i,j)/k1(i,j))^f ;
                end
            end
        else                                   % use linear interpolation
            kk2 = f*k2 + (1-f)*k1;             % 24 sept 04, k1, k2 switched
        end
        
        meank = (kk1+kk2)/2;                     % section k = mean of start and end values
        %  create element stiffness from section k
        
        vector = [1 1 1 1 1 1]*length;            % create H matrix
        H = diag(vector);
        H(4,2) = -(length^2)/2;
        H(5,1) = - H(4,2);
        Q = diag(vector)*length/2;    % create Q matrix
        Q(4,2) = -(length^3)/3;
        Q(5,1) = -Q(4,2);
        % form element stiffness matrix
        matrix = inv(inv(meank)*H + E*inv(meank)*Q);
        
    else
        %%  use Ansys results
        
        if feaflag == 0   % check for lack of fea model
            disp('error: no fea model but inadequate givenk coverage')
            break;
        end
        % create matrix of relative displacements
        for iload = 1:6  % loop on load cases
            relativedispl(1,iload) = sectiondispl(ielement+1,1,iload)...
                -sectiondispl(ielement,1,iload)-length*sectiondispl(ielement,5,iload);
            relativedispl(2,iload) = sectiondispl(ielement+1,2,iload)...
                -sectiondispl(ielement,2,iload)+length*sectiondispl(ielement,4,iload);
            relativedispl(3,iload) = sectiondispl(ielement+1,3,iload)...
                -sectiondispl(ielement,3,iload);
            relativedispl(4,iload) = sectiondispl(ielement+1,4,iload)...
                -sectiondispl(ielement,4,iload);
            relativedispl(5,iload) = sectiondispl(ielement+1,5,iload)...
                -sectiondispl(ielement,5,iload);
            relativedispl(6,iload) = sectiondispl(ielement+1,6,iload)...
                -sectiondispl(ielement,6,iload);
            
        end   % end of loop on load cases
        
        % create matrix of loads at end 2 of element
        
        force = zeros(6,6);
        for iload = 1:6    % loop on loads
            
            force(iload,iload) = 1;
            if iload == 1
                force(5,iload) = zbeamnode(nbeamnode) - zbeamnode(ielement+1);
            end
            if iload == 2
                force(4,iload) = -zbeamnode(nbeamnode) + zbeamnode(ielement+1);
            end
            
        end     % end of loop on loads
        
        % create stiffness matrix by dividing forces by displacements
        %  create flexibility matrix and make it symmetric
        
        method=2;  % flag is here in order to try out different matrix inverse methods.

        switch method
            case 1
                % original method
                matrix = force/relativedispl;
                flex = relativedispl/force;
            case 2
                % Nate's Method 1
                matrix = spde(relativedispl',force',false);
                flex = spde(force',relativedispl',false);
            case 3
                % Nate's Method 2
                matrix = spp(relativedispl',force',false);
                flex = spp(force',relativedispl',false);
        end
        
        flex = (flex+flex')/2;   % Q: Is this needed for methods 2 and 3?  guess it doesn't hurt
        
    end   % end of if for given k
    
    %  store matrices in one array for later filing
    
    i1 = (ielement-1)*8;
    storestiff(i1+1:i1+6,1:6) = matrix;
    matrix = (matrix + matrix')/2;   % make matrix symmetric
    storestiffsym(i1+1:i1+6,1:6) = matrix;
    storeflexsym(i1+1:i1+6,1:6) = flex;
    
    lambda = eig(matrix);
    eigstore(:,ielement) = lambda;  %extract eigenvalues and store columnwise
    mineig = min(lambda);
    if mineig <= 0.0                % if min e'value is <=0 then write warning and update singular array
        disp(['negative or zero eigenvalue in element #' num2str(ielement)])
        nsingular = nsingular + 1;
        singular(nsingular) = ielement;
        fprintf(fid6,'singular element: %3i\n',ielement);
    end
    
end  % end of loop on elements

%% save files

save stiffmatrixsym.txt storestiffsym -ascii  % save stiffness matrices
save flexmatrixsym.txt storeflexsym -ascii  % save stiffness matrices
save eigstore.txt eigstore -ascii  % save eigenvalues

storeData.storestiff = storestiff;
storeData.storestiffsym = storestiffsym;

singularData.nsingular = nsingular;
singularData.singular = singular;
singularData.fid6 = fid6;
