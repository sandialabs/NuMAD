function precomp = NuMAD2PreComp(data,matdb)

% First, filter the full set of material properties and stacks
% contained in MatDBsi.txt down to only the materials and stacks that are
% used by the actual NuMAD model, as described by *.nmd.  fprintf(fid,'lays filter
% results to the screen and the saves final set of materials and stacks to
% Mat and Comp variables.

% Prep materials information
matids=[];
compids=[];
fid=fopen('PrepMat.txt','wt');

for ii=2:length(data.active.list)
    name=data.active.list{ii};
    fprintf(fid,'****************\n');
    fprintf(fid,'%s\n',name);
    fprintf(fid,'        contains:\n');
    
    % find index to composite material
    for j=1:length(matdb)
        if strcmp(matdb(j).name,name)
            index1=j;
            break;
        end
    end
    compids=[compids index1];
    
    fprintf(fid,' LAYERS:\n');
    for j=1:length(matdb(index1).layer)
        name=matdb(index1).layer(j).layerName;
        fprintf(fid,'  %s\n',name);
        for k=1:length(matdb)
            if strcmp(matdb(k).name,name)
                index2=k;
                break;
            end
        end
        fprintf(fid,'      Material ID# %i\n',index2);
        matids=[matids index2];
    end
end

fclose(fid);

matkeepers=unique(matids);
compkeepers=unique(compids);

Mat=matdb(matkeepers);
Comp=matdb(compkeepers);

% check that all matl properties are in place
for i=1:length(Mat)
    % Applicable to the isotropic matls:
    if isempty(Mat(i).ey)
        Mat(i).ey=Mat(i).ex;
    end
    if isempty(Mat(i).gxy)
        Mat(i).gxy=Mat(i).ex/(2*(1+Mat(i).nuxy));
    end
    if isempty(Mat(i).prxy)
        Mat(i).prxy=Mat(i).nuxy;
    end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

N_sections=length(data.station);
Bl_length=data.station(end).LocationZ;

for i=1:N_sections
    precomp.station(i).span=data.station(i).LocationZ/Bl_length;
    precomp.station(i).lepos=data.station(i).Xoffset;
    precomp.station(i).chlen=data.station(i).Chord;
    precomp.station(i).twist=data.station(i).DegreesTwist;
    precomp.station(i).tetype=data.station(i).TEtype;
end

% Convert NuMAD Keypoints to PreComp airfoil coords and create shape files
for i=1:N_sections
    x=data.station(i).coords(:,1)';
    y=data.station(i).coords(:,2)';
    % convert to chord-normalized coordinated of the airfoil nodes
    % the first node, a leading-edge node, must be (0,0)
    % find le node by finding the smallest x coordinate (should be zero)
    [~,le(i)]=min(data.station(i).coords(:,1));
    x=[x(le(i):end) x(1:le(i)-1)];
    y=[y(le(i):end) y(1:le(i)-1)];
    x=x-x(1);
    y(1)=0;
    
    % if flatback:
    if strcmpi(data.station(i).TEtype,'flat')
        pointer=find(y==0);
        pointer=pointer(2);
        %remove trailing edge point and replace later with shear web
        x=x([1:pointer-1 pointer+1:end]);
        y=y([1:pointer-1 pointer+1:end]);
    end
    
    % The following are unnecessary if max of x is exactly 1?
    r=1/max(x);
    x=x*r;
    y=y*r;
    
    % check to see that points are not too close together
    % (could just resample the airfoil at this point using resampleAirfoil.m)
    newx=x(1);
    newy=y(1);
    for j=2:length(x)
        % ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % change this so flatbacks can have multiple points specified at
        % the TE position of x=1
% %         %dist=sqrt( (x(j)-newx(end))^2 + (y(j)-newy(end))^2 );
        dist=sqrt( (x(j)-newx(end))^2 );

        %        if dist>0.013  % this worked for SDM2010 paper
        if dist>0.015 % this worked for SDM2010 paper
            newx=[newx x(j)];
            newy=[newy y(j)];
        elseif x(j)==1 && newx(end)~=1 % ble
            newx=[newx(1:end-1) x(j)];
            newy=[newy(1:end-1) y(j)];
        elseif x(j)==1 && newx(end)==1 % ble
            newx=[newx x(j)];
            newy=[newy y(j)];
        end
        % ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    end
    x=newx;
    y=newy;
    
    % The following are unnecessary if max of x is exactly 1?
    r=1/max(x);
    x=x*r;
    y=y*r;
    
    % specify "blunt trailing edge" for flatback airfoils
    if strcmpi(data.station(i).TEtype,'flat')
        % ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % this is now unnecessary since the x=1 points are not removed by
        % "newx" variable distance sampling
% %         pointer=find(x>(1-0.000001)&x<(1+0.000001));
% %         if y(pointer)>0
% %             x(pointer+1)=1;
% %         elseif y(pointer)<0
% %             x(pointer-1)=1;
% %         end
        % ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    else       % make sure trailing edge is at 1,0 for other airfoils
        pointer=find(x>(1-0.000001)&x<(1+0.000001));
        y(pointer)=0;
    end
    
    % remove last point(s) to prevent errors in PreComp
    n_to_remove=1;
    x=x(1:end-n_to_remove);
    y=y(1:end-n_to_remove);
    
    % use this plot to check the airfoil shapes
    if 0
        figure(101)
        plot(x,y,'bo')
        axis equal
        title(['Section ' num2str(i)])
        grid on
        pause(0.5)
        %pause
    end
    
    precomp.station(i).x=x;
    precomp.station(i).y=y;
end


% Set up Materials

for i=1:length(Mat)
    precomp.material.id(i)=i;
    precomp.material.e1(i)=Mat(i).ex;
    precomp.material.e2(i)=Mat(i).ey;
    precomp.material.g12(i)=Mat(i).gxy;
    precomp.material.nu12(i)=Mat(i).prxy;
    precomp.material.dens(i)=Mat(i).dens;
    precomp.material.name(i)={Mat(i).name};
end

% Set up Layup information

for i=1:length(Comp)
    precomp.stack(i).name={Comp(i).name};
    if ~strcmpi(Comp(i).thicknessType,'Constant')
        warning('PreComp is not set up for anything other than ''Constant'' thickness type')
    end
    precomp.stack(i).n_layer=Comp(i).uniqueLayers;
    if ~strcmpi(Comp(i).symmetryType,'none')
        warning('PreComp is not set up for anything other than symmetry ''none''')
    end
    for j=1:precomp.stack(i).n_layer
        precomp.stack(i).layer(j).n=Comp(i).layer(j).quantity;
        precomp.stack(i).layer(j).thickness=Comp(i).layer(j).thicknessA;
        precomp.stack(i).layer(j).theta=Comp(i).layer(j).theta;
        precomp.stack(i).layer(j).n=Comp(i).layer(j).quantity;
                
        % determine the material id for each layer
        for k=1:length(precomp.material.id)
            if strcmp(Comp(i).layer(j).layerName,precomp.material.name{k})
                precomp.stack(i).layer(j).matid=k;
            end
        end
    end
end


% Set up Shear Web Layup information
% determine the composite id for each shear web
SW=data.shearweb;
if ~isempty(SW)
    for j=1:length(SW)
        for k=1:length(precomp.stack)
            if strcmp(SW(j).Material,precomp.stack(k).name)
                SW(j).stackid=k;
            end
        end
    end
    
    % Find the ends of shear webs
    % If the end of one web is the beginning of another, then it isn't
    % counted
    non_terminal_web_id=[];
    for j=1:length(SW)
        sw_end_dps=SW(j).Corner([3 4]);
        sw_end_at=SW(j).EndStation;
        for k=1:length(SW)
            if k~=j
                if sw_end_at==SW(k).BeginStation
                    % getting here means that the station where the jth web ends at is a beginning station for another web
                    % need to find out whether the web is not the same continuous web
                    if sw_end_dps(1)==SW(k).Corner(2) && sw_end_dps(2)==SW(k).Corner(1)
                        % Note that this shear web should be counted as a
                        % terminal web
                        non_terminal_web_id=[non_terminal_web_id j];
                    end
                end
            end
        end
    end
    
    % set up indices to the SW's that are terminal (meaning non-nonterminal SW's)
    terminal_web_id=[];
    for i=1:length(SW)
        if isempty(find(i==non_terminal_web_id))
            terminal_web_id=[terminal_web_id i];
        end
    end
else
    % what if SW doesn't exist (meaning there is no shear web)
    
end

% Set up station data
for i=1:N_sections
    % check to make sure there is a leading edge dp in the section.  If
    % not, send an error message
    if isempty(find(data.station(i).dp==0))
        error('Each station must contain a material D.P. at the leading edge in order for NuMAD2PreComp to work.');
    end
    hp_pointer=find(data.station(i).dp<=0);
    hp_pointer=hp_pointer(end:-1:1);
    lp_pointer=find(data.station(i).dp>=0);
    perChord_hp=data.station(i).dp(hp_pointer)*-1;
    perChord_lp=data.station(i).dp(lp_pointer);
    N_divs_hp=length(perChord_hp)-1;
    N_divs_lp=length(perChord_lp)-1;
    mat_hp=data.station(i).sm(N_divs_hp:-1:1);
    mat_lp=data.station(i).sm(N_divs_hp+1:end);

    isflatback=strcmpi(data.station(i).TEtype,'flat');
        
    
    if isflatback     % if flatback:
%         if ~strcmp(mat_lp(end),mat_hp(end))
%             warning(['Flatback isn''t defined of consistent material in section #' num2str(i) '.  Using material specification for bottom (HP) half of flatback.' ])
%         end
        mat_flat=mat_hp(end);
        mat_lp=mat_lp(1:end-1);
        mat_hp=mat_hp(1:end-1);
        perChord_lp=perChord_lp(1:end-1);
        perChord_hp=perChord_hp(1:end-1);
        N_divs_lp=length(perChord_lp)-1;
        N_divs_hp=length(perChord_hp)-1;
    end
    
    % Store in PreComp structure
    precomp.station(i).hp.segmentloc=perChord_hp;
    precomp.station(i).lp.segmentloc=perChord_lp;
    precomp.station(i).hp.n_segments=N_divs_hp;
    precomp.station(i).lp.n_segments=N_divs_lp;
    
    % determine stack id for the materials before saving them
    % LP Surface
    for j=1:length(mat_lp)
        for k=1:length(precomp.stack)
            if strcmp(mat_lp(j),precomp.stack(k).name)
                precomp.station(i).lp.stackid(j)=k;
            end
        end
    end
    % HP Surface
    for j=1:length(mat_hp)
        for k=1:length(precomp.stack)
            if strcmp(mat_hp(j),precomp.stack(k).name)
                precomp.station(i).hp.stackid(j)=k;
            end
        end
    end
    
    % if at tip of blade, assign same as next inboard station
    if i==N_sections
        precomp.station(i).hp.stackid=precomp.station(i-1).hp.stackid;
        precomp.station(i).lp.stackid=precomp.station(i-1).lp.stackid;
    end
    
    % Do Shear web conversion stuff
    tmp=[];
    tmp2=[];
    if ~isempty(SW) % if there is a shearweb
        % Sort out Shear Web placement and material information
        for j=1:length(SW) % count number of shear webs at this station and find index to information
            if SW(j).BeginStation==i
                tmp=[tmp j];
            end
        end
        
        % assign sw stack id's to appropriate stations
        for j=1:length(tmp)
            precomp.station(i).sw.stackid(j)=SW(tmp(j)).stackid;
        end
        
        % find location of shearweb, fraction of chord, & perpendicular to chord
        for j=1:length(tmp)
            temp=SW(tmp(j)).Corner+1;
            chord_hp=data.station(i).dp(temp(2))*-1;
            chord_lp=data.station(i).dp(temp(1));
            precomp.station(i).sw.chloc(j)=mean([chord_hp chord_lp]);
        end
        
        % identify any terminal shear webs
        for j=1:length(terminal_web_id)
            if SW(terminal_web_id(j)).BeginStation==i
                disp(' ')
                disp(['NuMAD station #' num2str(i) ' has the beginning of a terminal shear web with SW id #' num2str(terminal_web_id(j)) ])
                tmp2=[tmp2 terminal_web_id(j)];
            end
        end
        sw_end_id{i+1}=tmp2;
        
    end  % end SW stuff, add in terminal shear web information later after this loop on stations
    
    % Deal with flatbacks as shear webs
    if isflatback
        precomp.station(i).fb.flag=1;
        % find index to composite material        
        for k=1:length(precomp.stack)
            if strcmp(mat_flat,precomp.stack(k).name)
                precomp.station(i).fb.stackid=k;                
                break;  % stop looping and searching
            end
        end
        % ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % if at tip of blade, assign same as next inboard station
        if i==N_sections
            precomp.station(i).fb.stackid=precomp.station(i-1).fb.stackid;
        end
        % ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    else
        precomp.station(i).fb.flag=0;
        precomp.station(i).fb.stackid=0;
    end    
end

% Now calculate more information for terminal webs
if ~isempty(SW)
    for i=2:N_sections
        % figure out how many are present in the structure already
        try
            jj=length(precomp.station(i).sw.stackid);
        catch
            jj=0;
        end
        
        % assign sw stack id's to appropriate stations
        for j=1:length(sw_end_id{i})
            precomp.station(i).sw.stackid(jj+j)=SW( sw_end_id{i}(j) ).stackid;
        end
        
        % Calculate shear web locations for teminal webs
        tmp=sw_end_id{i};
        for j=1:length(tmp)
            temp=SW(tmp(j)).Corner+1;
            chord_hp=data.station(i).dp(temp(3))*-1;
            chord_lp=data.station(i).dp(temp(4));
            precomp.station(i).sw.chloc(j+jj)=mean([chord_hp chord_lp]);
        end
    end
end

% Make sure that webs are in order from LE to TE of blade
for i=1:length(precomp.station)
    try
        tmp=[precomp.station(i).sw.chloc ; precomp.station(i).sw.stackid]';
        tmp=sortrows(tmp);
        precomp.station(i).sw.chloc = tmp(:,1)';
        precomp.station(i).sw.stackid = tmp(:,2)';
    catch
        % do nothing
    end
end

% pass precomp data structure into the PreComp file creation functions

writePreCompShape(precomp,'shape.inp')

writePreCompMaterial(precomp,'materials.inp')

writePreCompLayup(precomp,'layup.inp')

writePreCompInput(precomp,'input.pci')

fclose all;  % if not already, close all files


end