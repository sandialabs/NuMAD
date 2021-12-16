function [fatigueDamage,plotFatigue]= calcFatigue(ansysBladeMaterials,IEC,Ltheta,LthetaPlus90,Mtheta,MthetaPlus90...
             ,binnedElements,materials,sections,elements,plateStrainsTheta,plateStrainsThetaPlus90,iSegment)
     %bladeMaterials read in to obtain fatigue exponents
     

        %Initialize damage variables, get section data, and get element strains
    
    numElem = numel(binnedElements);

    if numElem ==0 
        if iSegment ~=13
            error('"zwidth" is smaller than element size. Decrease element size or increase zwdith')
        else 
            fatigueDamage = zeros(1,10); % Only do this for the web because sometimes the 
            plotFatigue=[];              % web does not span the length of
            return                       % the blade
        end
    end
    
    analysis='LU';
    env='noEnv';
    temp='noTemp';
    mfg='basicFlaw';
    calca='test';
    loada='12_directions';
    SFs=round(setIEC5(analysis,env,temp,mfg,calca,loada),2);

%     analysis='LF';
%     loada='2_directions';
%     calcb ='measuredSlope';
%     loadb = 'Markov';
%     SFf=setIEC5(analysis,env,temp,mfg,calca,loada,calcb,loadb);
    
    SFf = SFs; 
    warning('SFf = SFs for BAR project since we are doing more than 2 load directions but less than 12')
    fatigueDamage=zeros(numElem,10); 

    plotFatigue=[];
    
    for i=1:numElem
        elemFD=0;
        elemFDlayer=0;
        elemFDmat=0;
        elemFDflap=0;
        elemFDedge=0;
        elNo=plateStrainsTheta(binnedElements(i),1); %or eplateStrainsThetaPlus90, equivalently
        sec_num=elements(elNo,6);
        nLayers=size(sections.layers{sections.secID==sec_num},1);
        
        
        coordSys='local';
        localFieldsTheta=extractFieldsThruThickness(plateStrainsTheta,sections,elements,ansysBladeMaterials,elNo,coordSys);
        localFieldsThetaPlus90=extractFieldsThruThickness(plateStrainsThetaPlus90,sections,elements,ansysBladeMaterials,elNo,coordSys);

        
        npts=numel(localFieldsTheta.(['element' num2str(elNo)]).x3);
        nptsPerLayer=2;
        layer=sort([1:npts/nptsPerLayer 1:npts/nptsPerLayer]);
        for ix3=1:npts
            matNumber=localFieldsTheta.(['element' num2str(elNo)]).matNumber(ix3);
            if ~isempty(ansysBladeMaterials(matNumber).m)
                %Determine mean and amplitude stress
                MthetaFactor=localFieldsTheta.(['element' num2str(elNo)]).sig11(ix3)/Mtheta;
                MthetaPlus90Factor=localFieldsThetaPlus90.(['element' num2str(elNo)]).sig11(ix3)/MthetaPlus90;

                LthetaWithFactor=Ltheta; %Initialize for every point
                LthetaPlus90WithFactor=LthetaPlus90; 

                LthetaWithFactor(1,:)=LthetaWithFactor(1,:)*MthetaFactor; %Means
                LthetaWithFactor(:,1)=abs(LthetaWithFactor(:,1)*MthetaFactor); %Amplitudes

                LthetaPlus90WithFactor(1,:)=LthetaPlus90WithFactor(1,:)*MthetaPlus90Factor; %Means
                LthetaPlus90WithFactor(:,1)=abs(LthetaPlus90WithFactor(:,1)*MthetaPlus90Factor); %Amplitudes
                
                m=ansysBladeMaterials(matNumber).m;
                XTEN=materials(matNumber).XTEN*SFs; %Unfactor the factored resistance
                XCMP=materials(matNumber).XCMP*SFs; %Unfactor the factored resistance

                % Determine the maximum number of cycles for failure based
                % on available fatigue failure criterion or from data

                %Calculate fatigue damage value, layer, and material for flap and edge cycles
                if isfinite(MthetaFactor)
                    switch IEC.fatigueCriterion
                        case 'Shifted Goodman'  %one case for now
                            layerFDtheta=shiftedGoodman(LthetaWithFactor,XTEN,XCMP,m,SFs,SFf);
                    end
                else
                    layerFDtheta=0;
                end
                
                if isfinite(MthetaPlus90Factor)
                    switch IEC.fatigueCriterion
                        case 'Shifted Goodman' %one case for now
                            layerFDthetaPlus90=shiftedGoodman(LthetaPlus90WithFactor,XTEN,XCMP,m,SFs,SFf);
                    end
                else
                    layerFDthetaPlus90=0;
                end                   
                layerFD=layerFDtheta+layerFDthetaPlus90;
                
                if layerFD>elemFD
                    elemFD=layerFD;
                    elemFDlayer=layer(ix3);
                    elemFDmat=matNumber;
                end
                if layerFDtheta>=elemFDflap
                    elemFDflap=layerFDtheta;
                    elemFDflapLayer=layer(ix3);
                    elemFDflapMat=matNumber;
                end
                if layerFDthetaPlus90>=elemFDedge 
                    elemFDedge=layerFDthetaPlus90;
                    elemFDedgeLayer=layer(ix3);
                    elemFDedgeMat=matNumber;
                end

                FDvalue=[elemFD elemFDflap elemFDedge];
                FDlayer=[elemFDlayer elemFDflapLayer elemFDedgeLayer];
                FDmat=[elemFDmat elemFDflapMat elemFDedgeMat];
                fatigueDamage(i,:)=[elNo FDvalue FDlayer FDmat];
           end
        end
%         plotFatigue=[plotFatigue;elements(binnedElements(i),1) FDvalue(1)];
        
    end
    %Output results in comma-delimited format
    %table(fatigue_damage(:,1),fatigue_damage(:,2),fatigue_damage(:,5),fatigue_damage(:,8))

    [~,imax]=max(fatigueDamage(:,2));
    fatigueDamage = fatigueDamage(imax,:);
        
%         
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%
%         for j=1:nLayers %Loop through layers
%             %Get material values
%             matNumber = sections.layers{sections.secID==sec_num}(j,3);
%             
%             if ~isempty(ansysBladeMaterials(matNumber).m)
%                 %Determine mean and amplitude stress
%                 MthetaFactor=stressesTheta(binnedElements(i),j+2)/Mtheta;
%                 MthetaPlus90Factor=stressesThetaPlus90(binnedElements(i),j+2)/MthetaPlus90;
% 
%                 LthetaWithFactor=Ltheta; %Initialize for every layer
%                 LthetaPlus90WithFactor=LthetaPlus90; 
% 
%                 LthetaWithFactor(1,:)=LthetaWithFactor(1,:)*MthetaFactor; %Means
%                 LthetaWithFactor(:,1)=abs(LthetaWithFactor(:,1)*MthetaFactor); %Amplitudes
% 
%                 LthetaPlus90WithFactor(1,:)=LthetaPlus90WithFactor(1,:)*MthetaPlus90Factor; %Means
%                 LthetaPlus90WithFactor(:,1)=abs(LthetaPlus90WithFactor(:,1)*MthetaPlus90Factor); %Amplitudes
%                 
%                 m=ansysBladeMaterials(matNumber).m;
%                 XTEN=materials(matNumber).XTEN*SFs; %Unfactor the factored resistance
%                 XCMP=materials(matNumber).XCMP*SFs; %Unfactor the factored resistance
%                 
% 
%                 % Determine the maximum number of cycles for failure based
%                 % on available fatigue failure criterion or from data
% 
%                 %Calculate fatigue damage value, layer, and material for flap and edge cycles
%                 if isfinite(MthetaFactor)
%                     switch IEC.fatigueCriterion
%                         case 'Shifted Goodman'  %one case for now
%                             layerFDtheta=shiftedGoodman(LthetaWithFactor,XTEN,XCMP,m,SFs,SFf);
%                     end
%                 else
%                     layerFDtheta=0;
%                 end
%                 
%                 if isfinite(MthetaPlus90Factor)
%                     switch IEC.fatigueCriterion
%                         case 'Shifted Goodman' %one case for now
%                             layerFDthetaPlus90=shiftedGoodman(LthetaPlus90WithFactor,XTEN,XCMP,m,SFs,SFf);
%                     end
%                 else
%                     layerFDthetaPlus90=0;
%                 end                   
% 
%                     
%                     
%                 layerFD=layerFDtheta+layerFDthetaPlus90;
% %                 if elNo==1950
% %                     keyboard
% %                 end
%                 
%                 if layerFD>elemFD
%                     elemFD=layerFD;
%                     elemFDlayer=j;
%                     elemFDmat=matNumber;
%                 end
%                 if layerFDtheta>=elemFDflap
%                     elemFDflap=layerFDtheta;
%                     elemFDflapLayer=j;
%                     elemFDflapMat=matNumber;
%                 end
%                 if layerFDthetaPlus90>=elemFDedge 
%                     elemFDedge=layerFDthetaPlus90;
%                     elemFDedgeLayer=j;
%                     elemFDedgeMat=matNumber;
%                 end
% 
%                 FDvalue=[elemFD elemFDflap elemFDedge];
%                 FDlayer=[elemFDlayer elemFDflapLayer elemFDedgeLayer];
%                 FDmat=[elemFDmat elemFDflapMat elemFDedgeMat];
%                 fatigueDamage(i,:)=[elNo FDvalue FDlayer FDmat];
%            end
%         end
% %         plotFatigue=[plotFatigue;elements(binnedElements(i),1) FDvalue(1)];
%         
%     end
%     %Output results in comma-delimited format
%     %table(fatigue_damage(:,1),fatigue_damage(:,2),fatigue_damage(:,5),fatigue_damage(:,8))
% 
%     [~,imax]=max(fatigueDamage(:,2));
%     fatigueDamage = fatigueDamage(imax,:);
end



