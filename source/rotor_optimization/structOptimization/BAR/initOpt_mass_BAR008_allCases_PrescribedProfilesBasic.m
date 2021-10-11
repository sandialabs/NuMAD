function [] = initOpt_mass_BAR008_allCases_noSparCapWidth(designNuMADfolder, bladeFile)
close all; 
cd 'C:\data\BAR008-148-mk0p2-s1-fiberglass';
%% ************************************************************************
% Define the starting blade file to read in design details
% *************************************************************************
designNuMADfolder = 'C:\data\BAR008-148-mk0p2-s1-fiberglass\NuMAD';

close all;
blade = BladeDef;
%blade.readYAML('BAR008_Ernesto_NewMaterials.yaml')
%blade.readYAML('BAR008_Ernesto_JP3.yaml')
%blade.readYAML('BAR008_SNL.yaml')
%blade.readYAML('BAR008_ply_thick.yaml')
%blade.readYAML('BAR008_SNL_Larc.yaml')
blade.readYAML('BAR00.yaml')


cd(designNuMADfolder)
blade.mesh=.45;
UB=1.1;
LB=0.8;

npts=6; %Number of control points along the span
bugDist=0.02;
    %%
sparCmp = find(contains({blade.components.name},'Spar'));
coreCmp = find(contains({blade.components.name},'S_filler'));  %%Core in skin (excluding webs)
webCmp = find(contains({blade.components.name},'Web'));  %%Core in skin (excluding webs)
webCoreCmpTemp = find(contains({blade.components(webCmp).name},'filler'));  %%Core in skin (excluding webs)
webSkinCmpTemp = find(contains({blade.components(webCmp).name},'skin'))
webCoreCmp=webCmp(webCoreCmpTemp);
webSkinCmp=webCmp(webSkinCmpTemp);

teReinfCmp = find(contains({blade.components.name},'TE_reinforcement'));
leReinfCmp = find(contains({blade.components.name},'LE_reinf'));
    
    for i = 2:11
        blade.components(i).name
        X=blade.components(i).cp(:,1);
        Y=blade.components(i).cp(:,2);
        
        r=linspace(X(1),X(end),npts);
        
        %Add an extra point so that the code updateBOM will work
        r=[r X(end)]; 
        r(end-1)=X(end)-bugDist;
        figure(1)
        subplot(4,4,i-1)
        plot(X,Y,'k')
        hold on
        title([int2str(i) ' ' blade.components(i).name])
        y=interp1(X,Y,r);
        blade.components(i).cp=[r' y'];     

    end

        %Web 
    for i=[webCoreCmp]
        blade.components(i).name
        X=blade.components(i).cp(:,1);
        Y=blade.components(i).cp(:,2);
        
        figure(1)
        subplot(4,4,i-1)
        plot(X,Y,'k')
        
        hold on
        rweb=[X(1) X(end)-bugDist X(end)];
        y=interp1(X,Y,rweb); %Half of the layers outboard

        blade.components(i).cp=[rweb' y'];
        title([int2str(i) ' ' blade.components(i).name])
        plot(blade.components(i).cp(:,1),blade.components(i).cp(:,2),'ok--')
         
    end 
    
    
    for i=[webSkinCmp]
        blade.components(i).name
        X=blade.components(i).cp(:,1);
        Y=blade.components(i).cp(:,2);
        
        figure(1)
        subplot(4,4,i-1)
        plot(X,Y,'k')
        
        hold on
        rweb=[X(1) X(1)+ 0.2 X(end)];
        y=interp1(X,Y,rweb); %Half of the layers outboard
        blade.components(i).cp=[rweb' y'];
        title([int2str(i) ' ' blade.components(i).name])
        plot(blade.components(i).cp(:,1),blade.components(i).cp(:,2),'ok--')
         
    end 

%Initialize x0 size
    ct = 0; 
    imatch=[2, 11; %Shell Skin
            12 14; %Aft web skin
            15 17]; %rear web skin
    ctMatch=1;
    for i=2:17
        if i~=imatch(ctMatch,2)
            ct = ct+ length(blade.components(i).cp)-1 ; %The minus one is due to the need to equate the last two points due to Numad bugs for Precomp
        else
            ctMatch=ctMatch+1;
        end
    end
    
    x0=zeros(ct,1);
    lb=zeros(ct,1);
    ub=zeros(ct,1);

%Initialize x0 and upper and lower bound values

ctMatch=1;
ct=0;
for i=2:17
    if i~=imatch(ctMatch,2)
        nvar = length(blade.components(i).cp);
        for j = 1:nvar
            if j< nvar
                ct = ct+1;
                x0(ct)=blade.components(i).cp(j,2);
                ub(ct)=UB*x0(ct);
                lb(ct)=LB*x0(ct);
            else
                x0(ct)=blade.components(i).cp(j-1,2);
                ub(ct)=UB*blade.components(i).cp(j-1,2);
                lb(ct)=LB*blade.components(i).cp(j,2);  
            end

        end
    else
        ctMatch=ctMatch+1;
    end
end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize x0 size
compIndex=[sparCmp,leReinfCmp,teReinfCmp,coreCmp];
    ct=numel(compIndex);
    ub=UB*ones(ct,1);
    lb=LB*ones(ct,1);
     x0=num2cell(zeros(ct,1));

%Plot upper and lower bounds
ct=1;
ctMatch=1;
for i=compIndex
    xplot=blade.components(i).cp(:,1);
    x0{ct}=blade.components(i).cp(:,2);
    figure(1)
    subplot(4,4,i-1)
    plot(xplot,blade.components(i).cp(:,2),'ok--')
    y=lb(ct)*blade.components(i).cp(:,2);
    y(y<0)=0;
    plot(xplot,y,'rd') 
    y=ub(ct)*blade.components(i).cp(:,2);
    plot(xplot,y,'bd') 
    ct=ct+1;
end

%ub(ub<1) = 1;
% set optimization dependencies
A=[ 0  0  -1  1  0;
    0  0  0  -1  1];
b = [0; 0];

% set genetic algorithm optimization iteration limits
% pop=36; gen=10;
pop=36; gen=2;
pop=numel(x0)*10;



%% ************************************************************************
% Call the optimization routine - first run or continuation of a run??
% *************************************************************************
useRestartFile = 0;     % flag to continue from a previous analysis or not
useParallel = true;     % flag to run the code in parallel or not

% finalize the geometry
blade.updateGeometry
blade.updateKeypoints
blade.updateBOM

if useRestartFile
    restartFile='layupDesignCandidates.txt';
else
    restartFile='';
end


structOpt_mass_BAR008_PrescribedProfilesBasic(blade,lb,ub,A,b,pop,gen,restartFile,useParallel,compIndex,x0)





