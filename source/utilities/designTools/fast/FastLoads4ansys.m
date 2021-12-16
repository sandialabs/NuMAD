function loadsTableSorted = FastLoads4ansys(output,fast_gage,IEC,varargin)

%% read in the FAST main files to determine span location of FAST gages
hm = pwd;
cd ..

fst=readFastMain([IEC.fstfn '.fst']);
ad=readFastAD(fst.ADFile(2:end-1));
bld=readFastBlade(fst.BldFile{1}(2:end-1));

cd(hm);
% set up the blade span array for FAST output gage locations
rGage = [0; ad.RNodes(fst.Out.BldGagNd)-fst.TurbConf.HubRad];
bladeLength = (fst.TurbConf.TipRad-fst.TurbConf.HubRad);

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

%% find the maximum moment distribution from the design load case set
for bb=1 % only using results from one blade currently
    for tt=1:length(fast_gage.theta{bb})
        gage_names = fast_gage.labels{bb}{tt}';
        for gg = 1:length(gage_names)
            for vv=1:length({output.(gage_names{gg})})
                normalFS = {'1p1' '1p2' '1p3' '1p4' '1p5' '2p1' '3p2' '3p3' '4p2' '5p1' '6p1' '6p3'};
                abnormalFS = {'2p2' '2p3' '6p2' '7p1'};
                % add the specified partial safety factor for loads
                % NOTE: not adding critical deflection partial safety factor
                if contains(output(vv).(gage_names{gg}).Name, normalFS)
                    FSloads = 1.35; % normal safety factor (IEC 61400-1)
                elseif contains(output(vv).(gage_names{gg}).Name, abnormalFS)
                    FSloads = 1.1; % abnormal safety factor (IEC 61400-1)
                else
                    error('design load case not recognized')
                end
                fastVar(vv,gg)=output(vv).(gage_names{gg}).data * FSloads;
            end
        end
        if contains(gage_names{1},'max','IgnoreCase',true)
            Mrb_i{bb}(tt,:) = max(fastVar,[],1);
        elseif contains(gage_names{1},'min','IgnoreCase',true)
            Mrb_i{bb}(tt,:) = min(fastVar,[],1);
        else
            error('Variable naming not compatible with maximum or minimum')
        end    
        thetaMomentRotation(tt,1) = fast_gage.theta{bb}{tt}(1);
    end
end

for bb = 1 % assign one maximum/minimum value based from all of the blades
    Mrb = Mrb_i{bb};
end

%% EMA added:
[FzAtrGage,MzAtrGage] = getMaxZForceMoment(output);
Fzb_FAST = interp1([rGage;bladeLength],[FzAtrGage 0],rBladeForce,'pchip');
Mzb_FAST = interp1([rGage;bladeLength],[MzAtrGage 0],rBladeMoment,'pchip');
%% END

%% Determine the Fxb and Fyb forces from the resultant moment maxima
for tt = 1:length(thetaMomentRotation)
    %Transform input resultant moment to the ANSYS coordinate system
    inverseRotationMatrix = [cosd(thetaMomentRotation(tt)) -1*sind(thetaMomentRotation(tt)); ...
        sind(thetaMomentRotation(tt)) cosd(thetaMomentRotation(tt))];
    
    Mb_FAST = inverseRotationMatrix * [Mrb(tt,:); zeros(size(Mrb(tt,:)))];
    
    Mxb_FAST = interp1([rGage; bladeLength],[Mb_FAST(1,:) 0],rBladeMoment,'pchip');
    Myb_FAST = interp1([rGage; bladeLength],[Mb_FAST(2,:) 0],rBladeMoment,'pchip');
    

    
    [Fxb_FAST_ec(tt,:),r_out]=getForcesFromMoments(rBladeMoment,Myb_FAST,bladeLength);
    [Fyb_FAST_ec(tt,:),r_out]=getForcesFromMoments(rBladeMoment,Mxb_FAST,bladeLength);
%     qFxb_FAST_ble(tt,:)=getForcesFromMoments_ble_dist(rBladeMoment,Myb_FAST,bladeLength);
%     qFyb_FAST_ble(tt,:)=getForcesFromMoments_ble_dist(rBladeMoment,Mxb_FAST,bladeLength);
    Fxb_FAST(tt,:)=getForcesFromMoments_ble_force(rBladeMoment,Myb_FAST,bladeLength);
    %% EMA original:
%     Fyb_FAST(tt,:)=getForcesFromMoments_ble_force(rBladeMoment,Mxb_FAST,bladeLength);
    %% changed to:
    %% Reverse the sign of X-moments going in to be consistent with right-hand coordinate system
    Fyb_FAST(tt,:)=getForcesFromMoments_ble_force(rBladeMoment,-Mxb_FAST,bladeLength);
    %% END

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Recheck Moment Dist
    rBlade = [rBladeMoment; bladeLength];
    zout=rBladeForce;



    A=zeros(numel(zout));
    for i=1:numel(zout)
        A(i,i:end)=(zout(i:end)-rBlade(i))';
    end


    
    
    % save the input loads for use in ansys Fatigue analysis
    loads_table{tt}.input.rGage = rGage';
    loads_table{tt}.input.Mrb = Mrb(tt,:)*1000; %Store Mrb for ANSYS fatigue. Convert from kNm to Nm

    loads_table{tt}.input.Mxb = Mb_FAST(1,:)*1000; %Convert from kNm to Nm
    loads_table{tt}.input.Myb = Mb_FAST(2,:)*1000; %Convert from kNm to Nm
    
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

%Sort the loads table for angles increasing. Starting from zero degrees
nDirections = length(loads_table);
theta = zeros(1,nDirections);
for iDirection = 1:nDirections %Start with second direction
   %Get the load direction angle
   if mod(iDirection,2)
       theta(iDirection) = fast_gage.theta{1}{iDirection}(1);
   else
       theta(iDirection) = fast_gage.theta{1}{iDirection}(1)+180;
   end
end
[theta,pointer]=sort(theta);

loadsTableSorted=cell(1,4);
for iDirection = 1:nDirections %Start with second direction
    loadsTableSorted{iDirection} =loads_table{pointer(iDirection)};
    loadsTableSorted{iDirection}.theta=theta(iDirection);
end

end



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
