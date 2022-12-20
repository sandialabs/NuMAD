function LF=getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber,avgFaceStress,app,mat,coreMatName)
    %This is a recursive function for LF. It appends to the list of LF for each
    %element that has a positive LF.  EC
    [m,~]=size(avgFaceStress);
    for i=1:m %For each element in the current section number (ansysSecNumber)
        elno=avgFaceStress(i,1); %Element number
        S11a=avgFaceStress(i,2); %N/m^2 
        S22a=avgFaceStress(i,3); %N/m^2
        S12a=avgFaceStress(i,7); %N/m^2
        %Ignoring other stresses for the time being.
        %if elno==305
            %disp('press pause')
            %pause(10)

            [lf,  phicr]=checkWrinkle([S11a;S22a;S12a],mat,app,coreMatName);
            if lf>=0
               LF=[LF;ansysSecNumber elno lf phicr]; %Append load factor to list of load factors
            end
        %end
    end



    function [lf, phicr]=checkWrinkle(S_alphaBeta,mat,app,coreMatName)
        % For a single finite element, given the average in-plane stresses
        % in a face-sheet of that element, compute the load factor for that
        % element. EC

        % lf    - scalar load factor for the element
        % phicr - an angle, degrees. The direction of wrinkling
        % S11a,S22a,S12a - respective average face sheet stress
        % mat - material object
        % app - blade data


        %Locate the face sheet
        cellMat={};
        for i=1:length(mat.layer) 
            cellMat = [cellMat; {mat.layer(i).layerName}];   %Create a cell array to use "find"    
        end

        kbalsa=find(strcmp(coreMatName,cellMat)==1);
        iLayer=1:(kbalsa-1); %Index of distinct materials in the top face  
        %ilayer=(kbalsa+1):numel(cellMat)); %Number of distinct materials in the bottom face
        matCore=app.matdb(find(strcmp(coreMatName,app.matlist)==1));

        if strcmp(matCore.type,'orthotropic')   
            Ec=matCore.ez; %Note that matCore.ez is invarient for in-plane rotations only 
        elseif strcmp(matCore.type,'isotropic')
            Ec=matCore.ex;
            Gc=Ec/(2*(1+matCore.nuxy));
        else
            errordlg(sprintf('Material  "%s"  not found in database.', matCore.type),'Error');
            error('Material type "%s" not found in database.', matCore.type);
        end
        %ilayer=(kbalsa+1):numel(cellMat)); %Number of distinct materials in the bottom face
        dangle=2; %Resolution for angle (degrees)
        N=180/dangle+1; %Number of angles to evaluate
        angle = 0; 

        invLF=zeros(N,1);
        %Apt=zeros(3,3);
        %Bpt=zeros(3,3);
        Dpt=zeros(3,3);

        %Find total height of facesheet
        h=0;
        for klay = 1:numel(iLayer)
            h = h+mat.layer(iLayer(klay)).thicknessA * mat.layer(iLayer(klay)).quantity;
        end

        for kang = 1:N
           if strcmp(matCore.type,'orthotropic')   
                Gc=1/(sind(angle)^2*(1/matCore.gyz)+cosd(angle)^2*(1/matCore.gxz)); %In-plate clockwise rotation of: angle 
           end
            %Asssuming all ply angles are zero
            R_sig=[cosd(angle)^2, sind(angle)^2, -2*sind(angle)*cosd(angle);  %In-plate clockwise rotation of: angle 
                   sind(angle)^2, cosd(angle)^2, 2*sind(angle)*cosd(angle)
                   sind(angle)*cosd(angle), -sind(angle)*cosd(angle), cosd(angle)^2-sind(angle)^2];

            z1=-h/2; %The coordinate location of the bottom of the face (with the origin in the midplane of the face sheet
            for klay = 1:numel(iLayer)
                z2=z1+mat.layer(iLayer(klay)).thicknessA * mat.layer(iLayer(klay)).quantity;

                matklay = app.matdb(find(strcmp(mat.layer(iLayer(klay)).layerName,app.matlist)==1));
                % Bulid Plane Stress reduced compliance matrix for each
                % layer
                %fprintf('z1 = %f z2 = %f mat = %f   %s\n',z1,z2,matListnumber,mat.layer(klay).layerName)

                %Entries common to either isotropic or orthotropic entries
                Se=zeros(3);
                Se(1,1)=1/matklay.ex;
                Se(1,3)=0; %Valid for orthotropic materials only
                Se(2,3)=0; %Valid for orthotropic materials only
                Se(3,1)=0; %Valid for orthotropic materials only
                Se(3,2)=0; %Valid for orthotropic materials only

                if strcmp(matklay.type,'orthotropic')   
                    Se(1,2)=-matklay.prxy/matklay.ex;
                    Se(2,1)=-matklay.prxy/matklay.ex;
                    Se(2,2)=1/matklay.ey;
                    Se(3,3)=1/matklay.gxy;
                elseif strcmp(matklay.type,'isotropic')
                    Se(1,2)=-matklay.nuxy/matklay.ex;
                    Se(2,1)=-matklay.nuxy/matklay.ex;
                    Se(2,2)=1/matklay.ex;
                    Se(3,3)=2*(1+matklay.nuxy)/matklay.ex;
                else
                    errordlg(sprintf('Material  "%s"  not found in database.', matklay.type),'Error');
                    error('Material type "%s" not found in database.', matklay.type);
                end

                %Apt=Apt+R_sig*inv(Se)*R_sig'*(z2-z1);

                %Bpt=Bpt+1/2*R_sig*inv(Se)*R_sig'*(z2^2-z1^2);
                Dpt=Dpt+1/3*R_sig*inv(Se)*R_sig'*(z2^3-z1^3);
                z1=z2;
            end 
            Pcr=-3/2*(2*Dpt(1,1)*Ec*Gc)^(1/3); %N/m
            Pphi=(R_sig(1,:)*S_alphaBeta)*h; %N/m

            invLF(kang)=Pphi/Pcr;
            angle=angle+dangle;
        end
        [invlf, phicr_i]=max(invLF);
        lf=1/invlf;

        phicr=(phicr_i-1)*dangle;

%         if lf>1e6
%             figure(2)
%             plot(angle, Pcr,'k')
%             xlabel('Angle, \phi [deg]')
%             ylabel('Load [N/m]')
%             hold on; 
% 
%             plot(angle, Pphi,'r')
%             plot(angle, lf*Pphi,'b')
%             legend('P_c_r','P_\phi',strcat(num2str(lf),'P_\phi'))
%             hold off;
%             figure(3)
%             plot(angle,invLF,'r')
%             hold on; 
%             plot(phicr,invlf,'*')
%             xlabel('Angle, \phi [deg]')
%             ylabel('1/\lambda [ ]')
%             hold off;
%             fprintf('\n %8.2g  %i | %8.2g %8.2g | %8.2g %8.2g %8.2g\n',lf,ansysSecNumber,Pcr(phicr_i),Pphi(phicr_i),S11a,S22a,S12a)
%             pause(10)
%         else
%             %fprintf('\n %8.2g  %i   %8.2g %8.2g   %8.2g %8.2g %8.2g\n',lf,ansysSecNumber,Pcr(phicr_i),Pphi(phicr_i),S11a,S22a,S12a)
%         end
    end
end