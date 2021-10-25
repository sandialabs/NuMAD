function [result]=extractFieldsThruThickness(fileName,sections,elements,ansysBladeMaterials,elementNumbers,coordSys)
  

    if ischar(fileName)
        %Used for reading element stresses
        pat='ELEM\s*ZCENT\s*EPS11\s*EPS22\s*EPS12\s*KAPA11\s*KAPA22\s*KAPA12\s*GAMMA13\s*GAMMA23';
        NCOLS=10;

        plateStrainsTheta = readANSYSElementTable(fileName,pat,NCOLS);
    else
        plateStrainsTheta=fileName;
    end

    numElem = numel(elementNumbers);

    nptsPerLayer=2;
    for i=1:numElem
        elNo=elementNumbers(i);
        elementStringName=['element' int2str(elNo)];
        sec_num=elements(elNo,6);
        layerData=sections.layers{sections.secID==sec_num}; %Important, flip the columns since code starts at the top layer.

        nLayers=size(layerData,1);
        k=find(plateStrainsTheta(:,1)==elNo);
        
        if isempty(k)
            error('Element not found in supplied file')
        end
        
        %Plate strains and curvatures
        eps=[plateStrainsTheta(k,3); plateStrainsTheta(k,4); plateStrainsTheta(k,5)];
        kappa=[plateStrainsTheta(k,6); plateStrainsTheta(k,7); plateStrainsTheta(k,8)];
        
        offset=sections.shellOffset{sections.secID==sec_num};
        h=sections.totalThick(sections.secID==sec_num); %Total thickness
        
        result.(elementStringName).x3=zeros(nptsPerLayer*nLayers,1);
        result.(elementStringName).eps11=zeros(nptsPerLayer*nLayers,1);
        result.(elementStringName).eps22=zeros(nptsPerLayer*nLayers,1);
        result.(elementStringName).eps33=zeros(nptsPerLayer*nLayers,1);
        result.(elementStringName).eps23=zeros(nptsPerLayer*nLayers,1);
        result.(elementStringName).eps13=zeros(nptsPerLayer*nLayers,1);
        result.(elementStringName).eps12=zeros(nptsPerLayer*nLayers,1);
        
        result.(elementStringName).sig11=deal(zeros(nptsPerLayer*nLayers,1));
        result.(elementStringName).sig22=deal(zeros(nptsPerLayer*nLayers,1));
        result.(elementStringName).sig12=deal(zeros(nptsPerLayer*nLayers,1));
        result.(elementStringName).matNumber=zeros(nptsPerLayer*nLayers,1);
        
        ct=1;
        if offset=='BOT'
            result.(elementStringName).x3(ct)=h;					   
        elseif offset=='MID'
            result.(elementStringName).x3(ct)=h/2;
        elseif offset=='TOP'
            result.(elementStringName).x3(ct)=0;
        end
        
        for j=nLayers:-1:1 %Loop through layers
            t=layerData(j,2);
            dx3=t/(nptsPerLayer-1);
            
            %Get material values
            matNumber = layerData(j,3);

            
            ansysBladeMaterials(matNumber).name
            angle=layerData(j,4); %Layup angle

            S0=zeros(6);
            %Orthotropic assumption
            S0(1,1)=1/ansysBladeMaterials(matNumber).ex;
            S0(1,2)=-ansysBladeMaterials(matNumber).prxy/ansysBladeMaterials(matNumber).ex;
            S0(1,3)=-ansysBladeMaterials(matNumber).prxz/ansysBladeMaterials(matNumber).ex;

            S0(2,1)=S0(1,2);
            S0(2,2)=1/ansysBladeMaterials(matNumber).ey;
            S0(2,3)=-ansysBladeMaterials(matNumber).pryz/ansysBladeMaterials(matNumber).ey;

            S0(3,1)=S0(1,3);
            S0(3,2)=S0(2,3);

            S0(4,4)=1/ansysBladeMaterials(matNumber).gyz;
            S0(5,5)=1/ansysBladeMaterials(matNumber).gxz;
            S0(6,6)=1/ansysBladeMaterials(matNumber).gxy;

            Se0=inPlaneMatrix(S0);
									  

            if strcmp(coordSys,'global')
                beta=[cosd(angle) -sind(angle) 0;    %In-plane clockwise rotation of: angle 
                     sind(angle) cosd(angle) 0;
                     0           0           1];

                Rsig=RsigmaMatrix(beta);
                Reps=inv(Rsig)';
                
                Rsige=inPlaneMatrix(Rsig);
                
                Q=Rsige*inv(Se0)*Rsige';
                
                %Out of plane strain
                S=Reps*S0*Reps';
                Set=outOfPlaneMatrix(S);
            elseif strcmp(coordSys,'local')
                Q=inv(Se0);
                
                %Out of plane strain;
                Set=outOfPlaneMatrix(S0);
            end

           for layerPoint=1:nptsPerLayer

                result.(elementStringName).x3(ct)
                epse=[eps(1)+result.(elementStringName).x3(ct)*kappa(1);
                     eps(2)+result.(elementStringName).x3(ct)*kappa(2);
                     eps(3)+result.(elementStringName).x3(ct)*kappa(3)];
                 
                if strcmp(coordSys,'local')
                    beta=[cosd(angle) sind(angle) 0;    %In-plane counter-clockwise rotation of: angle 
                         -sind(angle) cosd(angle) 0;
                         0           0           1];

                    Rsig=RsigmaMatrix(beta);
                    Reps=inv(Rsig)';
                    Repse=inPlaneMatrix(Reps);
                    epse=Repse*epse;
                end
                
               sige=Q*epse;
               result.(elementStringName).eps11(ct)=epse(1);

               result.(elementStringName).eps22(ct)=epse(2);
               result.(elementStringName).eps12(ct)=epse(3);
               
               result.(elementStringName).sig11(ct)=sige(1);
               result.(elementStringName).sig22(ct)=sige(2);
               result.(elementStringName).sig12(ct)=sige(3);
               
               %Out of plane strain
               epst=Set*Q*epse;
               
               result.(elementStringName).eps33(ct)=epst(1);
               result.(elementStringName).eps23(ct)=epst(2);
               result.(elementStringName).eps13(ct)=epst(3); 
               result.(elementStringName).matNumber(ct)=matNumber;
               %Set up next data point 
               result.(elementStringName).x3(ct+1)=result.(elementStringName).x3(ct)-dx3;
               ct=ct+1;
           end
           result.(elementStringName).x3(ct)=result.(elementStringName).x3(ct-1);
       
        end
        
        %Remove last point that was added at the end of the loop
        result.(elementStringName).x3=result.(elementStringName).x3(1:end-1);

    end
     



							
end  



