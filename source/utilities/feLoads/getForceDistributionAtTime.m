function loads_table = getForceDistributionAtTime(pointer,IEC,varargin)
% pointer.file - FAST filename
% pointer.time - time at which to get moments


%% read in the FAST main files to determine span location of FAST gages
hm = pwd;
cd ..

fst=readFastMain([IEC.fstfn '.fst']);
ad=readFastAD(fst.ADFile(2:end-1));
bld=readFastBlade(fst.BldFile{1}(2:end-1));

rGage = [0; ad.RNodes(fst.Out.BldGagNd)-fst.TurbConf.HubRad];
bladeLength = (fst.TurbConf.TipRad-fst.TurbConf.HubRad);

%% EMA original:
% out=loadFASTOutData(pointer.file);
%% changed to:
%% Get output data with forces and moments adjusted for structural twist
gageRot = interp1(bld.prop.BlFract.*bladeLength,bld.prop.StrcTwst,rGage(2:10));
out=loadFASTOutDataGageRot(pointer.file,gageRot);
%% END
[M] = momentsAtTime(out,pointer.time);
%% EMA added:
[F] = forcesAtTime(out,pointer.time);
%% END

cd(hm);
% set up the blade span array for FAST output gage locations

% set up the blade span array for ANSYS force stations
% NOTE: forces are at element centers, moments will be interpolated to
% Aerodyn element beginning point
if ~isempty(varargin)
    rBladeForce = varargin{1};
    if size(rBladeForce,2)>size(rBladeForce,1)
        rBladeForce = transpose(rBladeForce);
    end
    if rBladeForce(end) >= bladeLength || rBladeForce(1) == 0
        error('Blade span position is not possible')
    end
    rBladeMoment(1) = 0;
    elemSize(1) = rBladeForce(1)*2;    
    for ii = 2:length(rBladeForce)
        elemSize(ii) = (rBladeForce(ii) - (rBladeForce(ii-1)+elemSize(ii-1)/2) ) * 2;
        rBladeMoment(ii,1) = rBladeMoment(ii-1)+elemSize(ii);
    end
    if any(rBladeMoment >= rBladeForce) || any(elemSize <= 0)
        error('did not set up the blade span vector according to Aerodyn requirements')
    end
else % use original Aerodyn station locations
    rBladeForce = ad.RNodes-fst.TurbConf.HubRad;
    rBladeMoment = rBladeForce - ad.DRNodes./2;
end

prebend = interp1(bld.prop.BlFract.*bladeLength,bld.prop.PrecrvRef,rBladeForce);
presweep = interp1(bld.prop.BlFract.*bladeLength,bld.prop.PreswpRef,rBladeForce);

%% EMA added:
    Fzb_FAST = interp1([rGage; bladeLength],[F{3} 0],rBladeForce,'pchip');
    Mzb_FAST = interp1([rGage; bladeLength],[M{3} 0],rBladeMoment,'pchip');
%%

%% Determine the Fxb and Fyb forces from the resultant moments
    tt=1;
    Mxb_FAST = interp1([rGage; bladeLength],[M{1} 0],rBladeMoment,'pchip');
    Myb_FAST = interp1([rGage; bladeLength],[M{2} 0],rBladeMoment,'pchip');
    
%     m0=25e6;
%     Mxb_FAST = interp1([0; bladeLength],[0;0],rBladeMoment,'pchip');
%     Myb_FAST = interp1([0; bladeLength],[m0; m0],rBladeMoment,'pchip');
    
%     figure
%     hold on
%     plot(rBladeMoment,Mxb_FAST,'o-')
%     plot(rGage,M{1},'x--')
%     grid on
% %     Mx=Mi*sind(zThetaRotation); 
% %     My=Mi*cosd(zThetaRotation);
% %     diffz=diff(rBlade); 
% %     fx=getForcesFromMoments(diffz,My);
% %     fy=getForcesFromMoments(diffz,Mx);
    
    [Fxb_FAST_ec(tt,:),r_out]=getForcesFromMoments(rBladeMoment,Myb_FAST,bladeLength);
    [Fyb_FAST_ec(tt,:),r_out]=getForcesFromMoments(rBladeMoment,Mxb_FAST,bladeLength);
%     qFxb_FAST_ble(tt,:)=getForcesFromMoments_ble_dist(rBladeMoment,Myb_FAST,bladeLength);
%     qFyb_FAST_ble(tt,:)=getForcesFromMoments_ble_dist(rBladeMoment,Mxb_FAST,bladeLength);
    Fxb_FAST(tt,:)=getForcesFromMoments_ble_force(rBladeMoment,Myb_FAST,bladeLength);
    %% EMA original:
%     Fyb_FAST(tt,:)=getForcesFromMoments_ble_force(rBladeMoment,Mxb_FAST,bladeLength);
    %% changed to:
    %% Reverse sign of X-moments going in to be consistent with right-hand coordinate system
    Fyb_FAST(tt,:)=getForcesFromMoments_ble_force(rBladeMoment,-Mxb_FAST,bladeLength);
    %% END
%     disp('ad Forces')
%     Fxb_FAST(tt,:)'
%     Fyb_FAST(tt,:)'
%     pause(60)
%     
%     MxCheck = zeros(length(rBladeMoment),1);
%     MyCheck = MxCheck;
%     rCheck = rBladeMoment - 2.5;
%     for i1 = 1:length(rCheck)
%         for i2 = i1:length(rCheck)
%             dist = rBladeMoment(i2) - rCheck(i1);
%             MxCheck(i1) = MxCheck(i1) + Fyb_FAST(tt,i2)*dist;
%             MyCheck(i1) = MyCheck(i1) + Fxb_FAST(tt,i2)*dist;
%         end
%     end
%     disp('ad Moments')
%     MxCheck
%     MyCheck
%     pause(60)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Recheck Moment Dist
    rBlade = [rBladeMoment; bladeLength];
    zout=rBladeForce;

    %for i=1:numel(zout)
    %    zout(i)=mean([rBlade(i),rBlade(i+1)]);
    %end

    A=zeros(numel(zout));
    for i=1:numel(zout)
        A(i,i:end)=(zout(i:end)-rBlade(i))';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    figure(2000+tt)
    subplot(2,2,1)
    plot(rBladeMoment,-A*Fyb_FAST(tt,:)','ok-.',rGage,M{1},'xr-')
    grid on;
     hold on;
     legend('For ad2ANSYS','Original')
    title('Mx')
    subplot(2,2,2)
    hold on;
    diff(M{2})
    plot(rBladeMoment,A*Fxb_FAST(tt,:)','ok-.',rGage,M{2},'xr-')
    grid on;
     hold on;
     legend('For ad2ANSYS','Original')
    title('My')
    
    
    subplot(2,2,3)
    hold on;
    plot(r_out,Fyb_FAST(tt,:)','ok-.')
    grid on;
    title('Fy')

    subplot(2,2,4)
    hold on;
    plot(rBladeForce,Fxb_FAST(tt,:)','ok-.')
    grid on;
    title('Fx')
    
    hold on;
    
    % save the input loads for use in ansys Fatigue analysis
    loads_table{tt}.input.rGage = rGage';
    loads_table{tt}.input.Mxb = M{1}*1000; %Convert from kNm to Nm
    loads_table{tt}.input.Myb = M{2}*1000; %Convert from kNm to Nm
    %% EMA added:
    loads_table{tt}.input.Mzb = M{3}*1000;
    %% END
    
    % save the forces for use in layupDesign_ANSYSbuckling buckling analysis
    loads_table{tt}.rBlade = rBladeForce';
    loads_table{tt}.Fxb = Fxb_FAST(tt,:)*1000; %Convert from kN to N
    loads_table{tt}.Fyb = Fyb_FAST(tt,:)*1000; %Convert from kN to N
    %% EMA original:
%     % Centrifugal forces are not included - simulate ANSYS with centrifugal loads independently
%     loads_table{tt}.Fzb = zeros(size(rBladeForce'))*1000; %Convert from kN to N
    %% changed to:
    loads_table{tt}.Fzb = Fzb_FAST'*1000;
    %% END
    % Bending moments are used to generate the force-pressure distribution 
    % and not added to the ANSYS load set
    loads_table{tt}.Mxb = zeros(size(rBladeForce'))*1000; %Convert from kNm to Nm
    loads_table{tt}.Myb = zeros(size(rBladeForce'))*1000; %Convert from kNm to Nm
    %% EMA original:
%     % torsional moments should probably be added to the load set for
%     % ultimate strength failure calculations
%     loads_table{tt}.Mzb = zeros(size(rBladeForce'))*1000; %Convert from kNm  to Nm
    %% changed to:
    loads_table{tt}.Mzb = Mzb_FAST'*1000;
    %% END
    
    % save other information for transferring the loads to ANSYS
    loads_table{tt}.Alpha = zeros(size(rBladeForce'));
    loads_table{tt}.prebend = prebend';
    loads_table{tt}.presweep = presweep';    

end

function[M] = momentsAtTime(out,time)
    nGages = 10; %Assuming data is at 10 stations
    
    row=find(out.data(:,1)==time);
    %% EMA original:
    % M = {[zeros(1,nGages)];[zeros(1,nGages)]};
    %% changed to:
    M = {[zeros(1,nGages)];[zeros(1,nGages)];[zeros(1,nGages)]};
    %%END
    %Initialize
    bladeNo = 'b1';
    load='M';
    gageBaseName = ['Root' load];
    for i =1:nGages  %Assuming data is at 10 stations
        dir = 'x';
        labelName = [gageBaseName dir bladeNo];
        col = find(contains(out.list,labelName));
        M{1}(i)=out.data(row,col);
        
        dir = 'y';
        labelName = [gageBaseName dir bladeNo];
        col = find(contains(out.list,labelName));
        M{2}(i)=out.data(row,col);
        
        %% EMA added:
        dir = 'z';
        labelName = [gageBaseName dir bladeNo];
        col = find(contains(out.list,labelName));
        M{3}(i)=out.data(row,col);
        %% END
        
        gageBaseName = ['Spn' int2str(i) load 'L'];
    end
end

%% EMA added:

function[F] = forcesAtTime(out,time)
    nGages = 10; %Assuming data is at 10 stations
    
    row=find(out.data(:,1)==time); 
    F = {[zeros(1,nGages)];[zeros(1,nGages)];[zeros(1,nGages)]};
    %Initialize
    bladeNo = 'b1';
    load='F';
    gageBaseName = ['Root' load];
    for i =1:nGages  %Assuming data is at 10 stations
        dir = 'x';
        labelName = [gageBaseName dir bladeNo];
        col = find(contains(out.list,labelName));
        F{1}(i)=out.data(row,col);
        
        dir = 'y';
        labelName = [gageBaseName dir bladeNo];
        col = find(contains(out.list,labelName));
        F{2}(i)=out.data(row,col);
        
        dir = 'z';
        labelName = [gageBaseName dir bladeNo];
        col = find(contains(out.list,labelName));
        F{3}(i)=out.data(row,col);
        
        gageBaseName = ['Spn' int2str(i) load 'L'];
    end
end

%% END


%%
function py=getForcesFromMoments_ble_dist(rBlade,Mx,bladeLength)

py=zeros(size(Mx)); %distributed load (p) to be output
Nelem = length(rBlade);

for ii=Nelem:-1:1
    
    if ii == Nelem
        V2 = 0;
        M2 = 0;
        deltaZ = bladeLength-rBlade(ii);
    else
        V2 = Vnode(ii+1);
        M2 = Mx(ii+1);
        deltaZ = rBlade(ii+1)-rBlade(ii);
    end
    z1 = rBlade(ii);
    py(ii) = (Mx(ii)-M2+V2*deltaZ) / deltaZ^2
    Vnode(ii) = V2-py(ii)*deltaZ    
end
end


%%
function Fy=getForcesFromMoments_ble_force(rBlade,Mx,bladeLength)

Fy=zeros(size(Mx));
Nelem = length(rBlade);

for ii=Nelem:-1:1
    % assumes a force that is at the center of the element
    if ii == Nelem
        V2 = 0;
        M2 = 0;
        deltaZ = bladeLength-rBlade(ii);
    else
        V2 = Vnode(ii+1);
        M2 = Mx(ii+1);
        deltaZ = rBlade(ii+1)-rBlade(ii);
    end
    Fy(ii) = 2*(Mx(ii)-M2+V2*deltaZ) / deltaZ;
    Vnode(ii) = V2-Fy(ii);    
end
end

%%
function [Fy,zout]=getForcesFromMoments(rBlade,Mx,bladeLength)
rBlade = [rBlade; bladeLength];
zout=zeros(numel(rBlade)-1,1);

for i=1:numel(zout)
    zout(i)=mean([rBlade(i),rBlade(i+1)]);
end

A=zeros(numel(zout));
for i=1:numel(zout)
    A(i,i:end)=(zout(i:end)-rBlade(i))';
end
Fy=A\Mx;

end
