function designVar=layupDesign_ANSYSfatigue(ansysBladeMaterials,wt,rccdata,IEC,loadsTable,config)

    if any(contains(lower(config.ansys.analysisFlags.fatigue),'all'))
        nSegments=1;
    else
        nSegments=numel(config.ansys.analysisFlags.fatigue);
    end

    % Order of the segment names in segmentNamesReference
    % is very important. config.ansys.analysisFlags.fatigue can
    % be any order 
    segmentNamesReference=["HP_TE_FLAT","HP_TE_ReINF","HP_TE_PANEL", "HP_SPAR","HP_LE_PANEL","HP_LE","LP_LE",...
                    "LP_LE_PANEL","LP_SPAR","LP_TE_PANEL","LP_TE_REINF","LP_TE_FLAT"];
    nsegmentNamesReference=numel(segmentNamesReference);   

    markovSize=16;
    designVar={}; %Initialize
    Yr=IEC.designLife;

    %fst=readFastMain(['IEC_' IEC.fstfn '.fst']);
    %simtime=IEC.numSeeds*(fst.SimCtrl.TMax-IEC.delay); % simulated and rainflow counted time, seconds
    simtime=IEC.numSeeds*(IEC.SimTime-IEC.delay);
    nSpace=90/loadsTable{2}.theta; %assuming first position is zero degree and the next entry angular inrement throughout            
    nDirections = length(loadsTable);%Assuming all loadTables have the same number of directions
   %Get material, section (i.e. laminate), and strains for elements
    materials = readANSYSMatl('NuMAD/Materials.txt','NuMAD/Strengths.txt');
    sections = readANSYSSections('NuMAD/Sections.txt');
    elements = readANSYSElem(['NuMAD/Elements.txt']);
    
    for kTheta=1:nDirections/2  % Loop through half of the number of load direction 
                                 %since blade movements along single direction constitues 
                                 %two directions (e.g positive flap deflections and negative 
                                 %ones are two directions; both wich make
                                 %the flap cycles)
        loadsTableTheta=loadsTable{kTheta};
        theta=loadsTableTheta.theta;
        loadsTableThetaPlus90=loadsTable{kTheta+nSpace}; 
        
        nGage=numel(loadsTableTheta.input.rGage);
        gageNumber=(1:nGage)';
        [criticalElement,fatigueDamage,criticalLayerNo, criticalMatNo]=deal(zeros(nGage,1));
        criticalMat=cell(nGage,1);
        rGage=loadsTableTheta.input.rGage;
        
        MrTheta = loadsTableTheta.input.Mrb;
        MrThetaPlus90 = loadsTableThetaPlus90.input.Mrb;

        plotFatigue=[];
        for i=1:nSegments
            if any(contains(lower(config.ansys.analysisFlags.fatigue),'all'))
                title='All segments'; %Used for printing title for table.
                fileNameTheta=['plateStrains-all-' int2str(kTheta) '.txt'];
                fileNameThetaPlus90=['plateStrains-all-' int2str(kTheta+nSpace) '.txt'];
                iSegment=1;
            elseif ~ strcmpi(config.ansys.analysisFlags.fatigue(i),'webs')
                title=config.ansys.analysisFlags.fatigue(i); %Used for printing title for table.
                iSegment = find(strcmpi(segmentNamesReference,config.ansys.analysisFlags.fatigue(i))==1);
                fileNameTheta=['plateStrains-' int2str(iSegment) '-' int2str(kTheta) '.txt'];
                fileNameThetaPlus90=['plateStrains-' int2str(iSegment) '-' int2str(kTheta+nSpace) '.txt'];
            else
                title='Webs'; %Used for printing title for table.
                iSegment=nsegmentNamesReference+1;
                fileNameTheta=['plateStrains-' int2str(iSegment) '-' int2str(kTheta) '.txt'];
                fileNameThetaPlus90=['plateStrains-' int2str(iSegment) '-' int2str(kTheta+nSpace) '.txt']; 
            end     

            
            disp(fileNameTheta)
            disp(fileNameThetaPlus90)
            
            %Used for reading element stresses
            pat='ELEM\s*ZCENT\s*EPS11\s*EPS22\s*EPS12\s*KAPA11\s*KAPA22\s*KAPA12\s*GAMMA13\s*GAMMA23';
            NCOLS=10;

            plateStrainsTheta = readANSYSElementTable(fileNameTheta,pat,NCOLS);
            plateStrainsThetaPlus90 = readANSYSElementTable(fileNameThetaPlus90,pat,NCOLS);      
            
            
            for chSpan=1:nGage

                direction =int2str(theta);
                Ltheta=getMomentMarkov(rccdata,wt,Yr,simtime,markovSize,chSpan,direction);

                if theta+90 <180
                    direction =int2str(theta+90);
                else
                    direction =int2str(theta-90);
                end
                LthetaPlus90=getMomentMarkov(rccdata,wt,Yr,simtime,markovSize,chSpan,direction);


                Mtheta = interp1(rGage,MrTheta,rGage(chSpan));
                MthetaPlus90 = interp1(rGage,MrThetaPlus90,rGage(chSpan));

                zwidth=0.75; % [m] check fatigue for elements within a band of zwidth centered 
                           % at a blade gage location.
                z1=rGage(chSpan)-zwidth/2;
                z2=rGage(chSpan)+zwidth/2;

                binnedElements=intersect(find(plateStrainsTheta(:,2)<z2),find(plateStrainsTheta(:,2)>z1)); %Loop through elements in within zwidth centered at rGauge
    
                [fdData,plotFatigueChSpan] = calcFatigue(ansysBladeMaterials,IEC,Ltheta,LthetaPlus90,Mtheta,MthetaPlus90...
                    ,binnedElements,materials,sections,elements,plateStrainsTheta,plateStrainsThetaPlus90,iSegment);
                plotFatigue=[plotFatigue;plotFatigueChSpan];
                criticalElement(chSpan)=fdData(1);
                fatigueDamage(chSpan)=fdData(2);
                criticalLayerNo(chSpan)=fdData(5);
                criticalMatNo(chSpan)=fdData(8);
                criticalMat{chSpan}= ansysBladeMaterials(fdData(8)).name;

            end
        
    %         plotFatigueFileName=['plotFatigue-' int2str(kTheta)];
    %         writePlotFatigue(plotFatigueFileName,plotFatigue)
            fprintf('\n\n\n ************************ Segment No-%i: %s ************************\n', i,title)
         
            table(gageNumber,criticalElement,fatigueDamage,criticalLayerNo,criticalMatNo,criticalMat)
    %         designVar{end+1}=max(fatigueDamage);
            designVar{kTheta}.fatigueDamage(i,:)=fatigueDamage;
            designVar{kTheta}.criticalElement(i,:)=criticalElement;
            designVar{kTheta}.criticalLayerNo(i,:)=criticalLayerNo;
            designVar{kTheta}.criticalMatNo(i,:)=criticalMatNo;
            
        end
    end
   %delete stresses-*-*.txt;
end  

function writePlotFatigue(fname,plotFatigue)
    %Write fatigue damage for each element. ANSYS requires elements to be
    %sorted
    n=length(plotFatigue(:,1));
    fid=fopen([fname, '.txt'],'w+');
        fprintf(fid,'Element fatigueDamage\n');
        plotFatigue=sortrows(plotFatigue,1);
        for i=1:length(plotFatigue(:,1))
            fprintf(fid,'%8i  %6.5E\n',plotFatigue(i,1),plotFatigue(i,2));
        end
    fclose(fid);
    
    %Write plot commands
    fid=fopen([fname '.mac'],'w+');
    fprintf(fid,'/post1\n');
    fprintf(fid,'set,last\n');
    fprintf(fid,'plnsol,u,sum\n');

    fprintf(fid,'etab,test,u,X\n');

    fprintf(fid,'*get,max_e,elem,0,count\n');
    fprintf(fid,'*dim,d_res,array,max_e,3\n');
    fprintf(fid,'!Column 1 = Element Number\n');
    fprintf(fid,'!Column 2 = Where I''m putting result data\n');

    fprintf(fid,'*vget,d_res(1,1),elem,,elist\n');
    fprintf(fid,'*vfill,d_res(1,2),ramp,0,0\n');

    fprintf(fid,'*dim,d_results,array,%i,2   !Need to specify the same size array as the data being read in\n',n);
    fprintf(fid,'*vread,d_results(1,1),%s,txt,,jik,2,%i,,1\n',fname,n );
    fprintf(fid,'(F8.0,E13.5)\n');

    fprintf(fid,'*get,d_temp,parm,d_results,dim,x\n');

    fprintf(fid,'j=1\n');
    fprintf(fid,'i=1\n');

    fprintf(fid,'d_run=1\n');

    fprintf(fid,'*dowhile,d_run\n');

    fprintf(fid,'*if,d_res(i,1),EQ,d_results(j,1),THEN\n');
    fprintf(fid,'d_res(i,2)=d_results(j,2)\n');
    fprintf(fid,'j=j+1\n');
    fprintf(fid,'*endif\n');

    fprintf(fid,'*if,j,GT,d_temp,THEN\n');
    fprintf(fid,'d_run=0\n');
    fprintf(fid,'*endif\n');

    fprintf(fid,'i=i+1\n');
    fprintf(fid,'*enddo\n');
%     fprintf(fid,'j=1\n');
%     fprintf(fid,'*do,i,1,max_e\n');
%     fprintf(fid,'*if,j,LT,%i,THEN\n',n+1);
%     fprintf(fid,'*if,d_res(i,1),EQ,d_results(j,1),THEN\n');
%     fprintf(fid,'d_res(i,2)=d_results(j,2)\n');
%     fprintf(fid,'j=j+1\n');
%     fprintf(fid,'*endif\n');
%     fprintf(fid,'*endif\n');
%     fprintf(fid,'*enddo\n');


    fprintf(fid,'allsel,all\n');
    fprintf(fid,'*vput,d_res(1,2),elem,1,etab,test\n');
    fprintf(fid,'Pretab\n');

    fprintf(fid,'pletab,test\n');

    
end


