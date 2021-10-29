function BPEPost
%FUNCTIONNAME  One-line function description
% **********************************************************************
% *           Part of the SNL Wind Turbine Analysis Toolbox            *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% **********************************************************************
%   output = functionName(input)
%   Detailed function description describing inputs, outputs, and other necessary information.
%
%      input  = description of inputs
%  
%      output = description of outputs
%

%===== CREDITS & CHANGELOG ================================================
% 2011.06.09  brr: initial release
% yyyy.mm.dd  initials: description

load('numaddata')
fid2 = fopen('bpe.txt','r'); % opens basic input file
if fid2==-1, error('BPE input file "bpe.txt" not found'), end
nbeamnode = str2num(scandm2(fid2)); % no of beam nodes desired
string = scandm2(fid2); % spanwise coords of beam nodes
[zbeamnode nbeamnode] = sscanf(string,'%f,',nbeamnode); % puts node values into vector
string = scandm2(fid2); % precurve offsets
fclose(fid2);
[precurve number] = sscanf(string,'%f,',nbeamnode); % puts precurve offsets into vector

N=nbeamnode-1;
a=txt2mat('kout.txt');

% kout.txt contents:
% row 1 or (N-1)*9+1: ielement, massperlen, cofpyoffset, Ixxperlen, Iyyperlen, Izzperlen, Ixyperlen
% row 2 or (N-1)*9+2: cgxoffset, cgyoffset,    maxchord,  chordrot, LExoffset, LEyoffset
% row 3 thru 8 or (N-1)*9+3 thru +8: 6x6 section stiffness matrix
% row 9 or (N-1)*9+9: blank

if strmatch(data.BladeRotation,'ccw');rotdir=1;else;rotdir=-1;end

x_el_off=zeros(N,1);
y_el_off=zeros(N,1);
x_sh_off=zeros(N,1);
y_sh_off=zeros(N,1);
rot=zeros(N,1);
EI_flap=zeros(N,1);
EI_edge=zeros(N,1);
GJ=zeros(N,1);
Alphaa=zeros(N,1);
sh_flap=zeros(N,1);
sh_edge=zeros(N,1);
locationz=zbeamnode;
zbeamnode_mid1=mean([locationz(1:end-1) locationz(2:end)],2);
zbeamnode_mid=[0 ;zbeamnode_mid1 ;locationz(end)];
rot_aero=interp1(locationz,tw_aero,zbeamnode_mid1)/(rotdir*-57.3);

for j=1:N
    % gather mass and inertia information
    rows=9*(j-1)+1:9*(j-1)+2;
    M=a(rows,:);
    massperlen(j)=M(1,2);
    Ixx=M(1,4);
    Iyy=M(1,5);
    Ixy=M(1,7);
    cgi=[M(2,1) M(2,2)]';
    
    % rotation cg's to chord twist axis
    R=[cos(rot_aero(j)) sin(rot_aero(j));
        -sin(rot_aero(j)) cos(rot_aero(j))];
    cg(:,j)=R*cgi;
    
    % Rotate global inertias (blade pitch axis) to blade twist axis
    Ixx_rot=(Ixx+Iyy)/2+(Ixx-Iyy)/2*cos(2*rot_aero(j))-Ixy*sin(2*rot_aero(j));
    Iyy_rot=(Ixx+Iyy)/2-(Ixx-Iyy)/2*cos(2*rot_aero(j))+Ixy*sin(2*rot_aero(j));
    % Translate rotated inertias to axis at center of mass
    Ixxperlen_cg(j)=Ixx_rot-massperlen(j)*cg(2)^2;
    Iyyperlen_cg(j)=Iyy_rot-massperlen(j)*cg(1)^2;
    
    % gather stiffness information
    rows=9*(j-1)+3:9*(j-1)+8;
    A=a(rows,1:6);
    
    % x & y offsets to elastic axes from global c.s. (pitch axis)
    x_el_off(j)=-1*A(5,3)/A(3,3);
    y_el_off(j)=A(4,3)/A(3,3);
    x_sh_off(j)=A(2,6)/A(2,2);
    y_sh_off(j)=-1*A(1,6)/A(1,1);
    
    % x & y offsets to elastic axes from section c.s. (centered on pitch axis and rotated)
    el=[x_el_off(j) y_el_off(j)]';
    sh=[x_sh_off(j) y_sh_off(j)]';
    el=R*el;
    sh=R*sh;
    x_el_off_sect(j)=el(1);
    y_el_off_sect(j)=el(2);
    x_sh_off_sect(j)=sh(1);
    y_sh_off_sect(j)=sh(2);
    
    % k matrix after partial translation to elastic axis
    B=A;
    B(:,4)=A(:,4)-y_el_off(j)*A(:,3);
    B(:,5)=A(:,5)+x_el_off(j)*A(:,3);
    B(:,6)=A(:,6)+y_el_off(j)*A(:,1)-x_el_off(j)*A(:,2);
    
    % k matrix after partial translation to shear center
    D=A;
    D(:,4)=A(:,4)-y_sh_off(j)*A(:,3);
    D(:,5)=A(:,5)+x_sh_off(j)*A(:,3);
    D(:,6)=A(:,6)+y_sh_off(j)*A(:,1)-x_sh_off(j)*A(:,2);
    
    % k matrix wrt initial coordinates passing through elastic axis
    C=B;
    C(4,:)=B(4,:)-y_el_off(j)*B(3,:);
    C(5,:)=B(5,:)+x_el_off(j)*B(3,:);
    C(6,:)=B(6,:)+y_el_off(j)*B(1,:)-x_el_off(j)*B(2,:);
    
    EA2(j)=C(3,3);
    
    % k matrix wrt initial coordinates passing through shear axis
    E=D;
    E(4,:)=D(4,:)-y_sh_off(j)*D(3,:);
    E(5,:)=D(5,:)+x_sh_off(j)*D(3,:);
    E(6,:)=D(6,:)+y_sh_off(j)*D(1,:)-x_sh_off(j)*D(2,:);
    
    % rotation from initial to principle axes, radians
    if(abs((C(4,4)-C(5,5))/C(4,5))>0.01)
        rot_i(j)=0.5*atan(2*C(4,5)/(C(4,4)-C(5,5)));
    else
        rot_i(j)=0;
    end
    
    rot=rot_aero;  % use aerodynamic twist as axis for rotation transformation
    %rot=rot_i; % use principle axis for rotation transformation
    
    % rotation transformation to principle axes
    matrot=[cos(rot(j))  sin(rot(j)) 0 0            0             0;
        -sin(rot(j))     cos(rot(j)) 0 0            0             0;
        0                   0        1 0            0             0;
        0                   0        0 cos(rot(j))  sin(rot(j))   0;
        0                   0        0 -sin(rot(j)) cos(rot(j))   0;
        0                   0        0 0            0             1];
    
    % Section stiffness matrix, k, (wrt principal coordinates passing through elastic axis)
    sectstiffmatrix=matrot*C*matrot';
    % k matrix after translation to shear center and rotation to principal axes
    F=matrot*E*matrot';
    
    % inverse of k wrt principal coords passing through elastic axiis (i.e. flexibility matrix, S)
    sectflexelmatrix=inv(sectstiffmatrix);
    sectflexshmatrix=inv(F);
    
    % Flapwise pseudo EI [=1/S(5,5)], elastic axes
    if sectflexelmatrix(5,5)>0
        EI_flap(j)=1/sectflexelmatrix(5,5);
    else
        EI_flap(j)=sectstiffmatrix(5,5);
    end
    
    % Edgewise pseudo EI [=1/S(4,4)], elastic axes
    if sectflexelmatrix(4,4)>0
        EI_edge(j)=1/sectflexelmatrix(4,4);
    else
        EI_edge(j)=sectstiffmatrix(4,4);
    end
    
    % Torsional pseudo GJ [=1/S(6,6)], about shear center
    if sectflexshmatrix(6,6)>0
        GJ(j)=1/sectflexshmatrix(6,6);
    else
        GJ(j)=F(6,6);
    end
    
    % Flap-twist coupling coeff [=-S(5,6)/sqrt(S(5,5)*S(6,6))], shear center
    if sectflexshmatrix(6,6)*sectflexshmatrix(5,5)>0
        Alphaa(j)=-sectflexshmatrix(5,6)/sqrt(sectflexshmatrix(6,6)*sectflexshmatrix(5,5));
    else
        Alphaa(j)=F(5,6)/sqrt(F(6,6)*F(5,5));
    end
    
    % Flapwise shear stiffness [=k(1,1)], shear center
    sh_flap(j)=F(1,1);
    
    % Edgewise shear stiffness [=k(2,2)], shear center
    sh_edge(j)=F(2,2);
    
    % Flapwise pseudo EI [=1/S(5,5)], elastic axes
    if sectflexelmatrix(3,3)>0
        EA(j)=1/sectflexelmatrix(3,3);
    else
        EA(j)=sectstiffmatrix(3,3);
    end
    
end

xcg=cg(1,:)';
ycg=cg(2,:)';

% Extrapolation to end points, 0 and 1
tmp=log(EI_flap);
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
EI_flap=exp(tmp2);

tmp=log(EI_edge);
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
EI_edge=exp(tmp2);

tmp=log(GJ);
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
GJ=exp(tmp2);

tmp=Alphaa;
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
Alphaa=tmp2;

tmp=log(sh_flap);
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
sh_flap=exp(tmp2);

tmp=log(sh_edge);
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
sh_edge=exp(tmp2);

tmp=-1*rot_i';
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
rot_i=tmp2;

tmp=x_el_off_sect;
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
x_el_off_sect=tmp2;

tmp=y_el_off_sect;
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
y_el_off_sect=tmp2;

tmp=x_sh_off_sect;
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
x_sh_off_sect=tmp2;

tmp=y_sh_off_sect;
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
y_sh_off_sect=tmp2;

tmp=log(EA');
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
EA=exp(tmp2);

tmp=log(massperlen');
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
massperlen=exp(tmp2);

tmp=log(Ixxperlen_cg');
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
Ixxperlen_cg=exp(tmp2);

tmp=log(Iyyperlen_cg');
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
Iyyperlen_cg=exp(tmp2);

tmp=xcg;
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
xcg=tmp2;

tmp=ycg;
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
ycg=tmp2;

tmp=rot_aero;
tmp2=interp1(zbeamnode_mid1,tmp,zbeamnode_mid,'linear','extrap');
rot_aero=tmp2;


% Plotting
figure('Name','BPE Outputs Inspection')
subplot(3,4,1)
plot(zbeamnode_mid,EI_flap,'-o')
ylabel('EI\_flap [Nm^2]')
subplot(3,4,2)
plot(zbeamnode_mid,EI_edge,'-o')
ylabel('EI\_edge [Nm^2]')
subplot(3,4,3)
plot(zbeamnode_mid,GJ,'-o')
ylabel('GJ [Nm^2]')
subplot(3,4,4)
plot(zbeamnode_mid,Alphaa,'-o')
ylabel('Alphaa')
subplot(3,4,5)
plot(zbeamnode_mid,sh_flap,'-o')
ylabel('sh\_flap [N]')
subplot(3,4,6)
plot(zbeamnode_mid,sh_edge,'-o')
ylabel('sh\_edge [N]')
subplot(3,4,7)
plot(zbeamnode_mid,rot_i*57.3,'-o')
ylabel('rotation [^o]')
subplot(3,4,8)
plot(zbeamnode_mid,x_el_off_sect,'-o')
ylabel('x\_el\_off [m]')
subplot(3,4,9)
plot(zbeamnode_mid,y_el_off_sect,'-o')
ylabel('y\_el\_off [m]')
subplot(3,4,10)
plot(zbeamnode_mid,x_sh_off_sect,'-o')
ylabel('x\_sh\_off [m]')
subplot(3,4,11)
plot(zbeamnode_mid,y_sh_off_sect,'-o')
ylabel('y\_sh\_off [m]')
subplot(3,4,12)
plot(zbeamnode_mid,EA,'-o')
ylabel('EA [Nm]')

save BPE_SectionData zbeamnode_mid EI_flap EI_edge GJ Alphaa sh_flap sh_edge rot_i x_el_off_sect  y_el_off_sect x_sh_off_sect y_sh_off_sect EA massperlen Ixxperlen_cg Iyyperlen_cg xcg ycg locationz rot_aero rotdir

end

