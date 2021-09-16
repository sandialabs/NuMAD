function  Eb31out(nodeData,sectionData,flagData,storeData,massData,xyData,precurve)
%  organizes output to wtprep and for xls spreadsheet
%  15 nov  2002 mass per unit length added to xlsout
%  14 march 03, precurve offsets added
%  30 june2003 additional output file for givenk and given  inertias
%  latest 2 july 2003
%  20 nov 2003, v27, for use with Legacy Numad
%  26 jan 2004, v28, for use with Legacy Numad
%  30 march 2004, v29, mods for local coordinate option
%  13 april 2004, v30, outer profile included
%  22 july 2005, v31, outer profile removed

feaflag = flagData.feaflag;
convertflag = flagData.convertFlag;

storestiffsym = storeData.storestiffsym;
storesectstiff = storeData.storesectstiff;

nbeamnode = nodeData.nbeamnode;
zbeamnode = nodeData.zbeamnode;

nodemass = massData.nodemass;
cgx = massData.cgx;
cgy = massData.cgy;
cgz = massData.cgz;
massmtx2 = massData.massmtx2;
massmty2 = massData.massmty2;
massmtz2 = massData.massmtz2;
massmtxy = massData.massmtxy;
massmtxz = massData.massmtxz;
massmtyz = massData.massmtyz;
Ixx = massData.Ixx;
Iyy = massData.Iyy;
Ixy = massData.Ixy;

xlocate = xyData.xlocate;
ylocate = xyData.ylocate;
xrotate = xyData.xrotate;
yrotate = xyData.yrotate;

sectionnodecount1 = sectionData.sectionnodecount1;
sectionnodecoord1 = sectionData.sectionnodecoord1;
sectionnodecount2 = sectionData.sectionnodecount2;
sectionnodecoord2 = sectionData.sectionnodecoord2;
sectioncofpy2 = sectionData.sectioncofpy2;
sectioncofpx = sectionData.sectioncofpx;
sectioncofpy = sectionData.sectioncofpy;
sectionLE = sectionData.sectionLE;
sectionmaxchord = sectionData.sectionmaxchord;
sectionchordrot = sectionData.sectionchordrot;

%
[n1sectionnodecoord2,n2sectionnodecoord2,n3sectionnodecoord2]=size(sectionnodecoord2);

%  output for wtprep
%
fid3 = fopen('.\admout.txt','w');      % opens file with name of file for blade output info
%
presweep = zeros(nbeamnode,1);
temp1 = zeros(3,1);
temp2 = zeros(3,1);
temp = zeros(3,1);

if convertflag == 1
    for i = 1:nbeamnode                      % loop on nodes / elements for conversion to local axes
        cgx(i) = cgx(i) - xlocate(i);
        cgy(i) = cgy(i) - ylocate(i);
        sectioncofpx(i) = sectioncofpx(i) - xlocate(i);
        sectioncofpy(i) = sectioncofpy(i) - ylocate(i);
        sectionLE(i,1) = sectionLE(i,1) - xlocate(i);
        sectionLE(i,2) = sectionLE(i,2) - ylocate(i);
        precurve(i) = xlocate(i);
        presweep(i) = ylocate(i);
        
        R = zeros(3,3);
        R(1,1) = cos(yrotate(i));         % transformation matrix.  Note, R is for following element.
        R(1,2) = 0.0;
        R(1,3) = -sin(yrotate(i));
        R(2,1) = sin(xrotate(i))*sin(yrotate(i));
        R(2,2) = cos(xrotate(i));
        R(2,3) = sin(xrotate(i))*cos(yrotate(i));
        R(3,1) = cos(xrotate(i))*sin(yrotate(i));
        R(3,2) =-sin(xrotate(i));
        R(3,3) = cos(xrotate(i))*cos(yrotate(i));
        
        temp = R*[cgx(i);cgy(i);cgz(i)];    % transform cg offsets
        cgx(i) = temp(1);
        cgy(i) = temp(2);
        cgz(i) = temp(3);
        
        temp = R*[sectioncofpx(i);sectioncofpy(i);0];  % transform c of p offsets
        sectioncofpx(i) = temp(1);
        sectioncofpy(i) = temp(2);
        
        temp = R*[sectionLE(i,1);sectionLE(i,2);0];    % transform LE offsets
        sectionLE(i,1) = temp(1);
        sectionLE(i,2) = temp(2);
        
        temp1(1) = sectionmaxchord(i)*sin(sectionchordrot(i));
        temp1(2) = sectionmaxchord(i)*cos(sectionchordrot(i));
        temp1(3) = 0;
        temp = R*temp1;
        dx = temp(1);
        dy = temp(2);
        sectionmaxchord(i) = sqrt(dx^2 + dy^2);
        sectionchordrot(i) = atan(dx/dy);
        
    end                                      % end of loop on nodes / elements
else
    xlocate = zeros(nbeamnode,1);
    ylocate = zeros(nbeamnode,1);
    
end                                        % end of if conversion to local axes

for ielement = 1:nbeamnode-1;                    % loop on beam elements
    matrix = zeros(6,6);
    
    if feaflag==1                    % skip this section if fea model present
        
        mass = nodemass(ielement);
        maxchord = sectionmaxchord(ielement);
        chordrot = sectionchordrot(ielement);
        cofpxoffset = sectioncofpx(ielement);
        cofpyoffset = sectioncofpy(ielement);
        precurveoffset = precurve(ielement);
        presweepoffset = presweep(ielement);
        cgxoffset = cgx(ielement);
        cgyoffset = cgy(ielement);
        cgzoffset = cgz(ielement);
        
    else                             % do this if no fea model
        
        mass = nodemass(ielement);
        maxchord = sectionmaxchord(ielement);
        chordrot = sectionchordrot(ielement);
        cofpxoffset = sectioncofpx(ielement);
        cofpyoffset = sectioncofpy(ielement);
        precurveoffset = precurve(ielement);
        presweepoffset = presweep(ielement);
        cgxoffset = cgx(ielement);
        cgyoffset = cgy(ielement);
        cgzoffset = cgz(ielement);
    end
    %
    fprintf(fid3,'%9.3f %9.3f %9.2f %9.3f %9.3f %9.3f %9.3f !root_dist, chord, chord_rot, c_of_p offsets, precurve, presweep \n',...
        zbeamnode(ielement),maxchord,chordrot*57.3,cofpxoffset,cofpyoffset,precurveoffset,presweepoffset);
    fprintf(fid3,'%8.2f %8.3f %8.3f %8.3f  !mass, cgoffsets \n',...
        mass,cgxoffset,cgyoffset,cgzoffset);
    fprintf(fid3,'%11.3f %11.3f %11.3f %11.3f %11.3f %11.3f  !inertia tensor\n',...
        massmtx2(ielement),massmty2(ielement),massmtz2(ielement),massmtxy(ielement),massmtxz(ielement),massmtyz(ielement));
    matrix = storestiffsym((ielement-1)*8+1:(ielement-1)*8+6,1:6);
    fprintf(fid3,'%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n',matrix);
    
end;                                         % end of loop on beam elements

% print one line for end of model
if feaflag==1                    % skip this section if no fea model
    maxchord = sectionmaxchord(nbeamnode);
    chordrot = sectionchordrot(nbeamnode);
    cofpxoffset = sectioncofpx(nbeamnode);
    cofpyoffset = sectioncofpy(nbeamnode);
    precurveoffset = precurve(nbeamnode);
else                             % do this if no fea model
    maxchord = sectionmaxchord(nbeamnode-1);
    chordrot = sectionchordrot(nbeamnode-1);
    cofpxoffset = sectioncofpx(nbeamnode-1);
    cofpyoffset = sectioncofpy(nbeamnode-1);
    precurveoffset = precurve(nbeamnode-1);
    presweepoffset = presweep(nbeamnode-1);
end

fprintf(fid3,'%9.3f %9.3f %9.2f %9.3f %9.3f %9.3f %9.3f !root_dist, chord, chord_rot, c_of_p offsets, precurve, presweep \n',...
    zbeamnode(nbeamnode),maxchord,chordrot*57.3,cofpxoffset,cofpyoffset,precurveoffset,presweepoffset);

status=fclose(fid3);

%   output for spreadsheet  and givenk

fid4 = fopen('.\xlsout.txt','w');      % opens file with name of file for blade output info for spreadsheet
fid5 = fopen('.\kout.txt','w');      % opens file with name of file for blade givenk and given inertias

if feaflag == 1;                         % skip next section if no fea model
    
    %                                   % combine sections # 1 and 2
    sectionnodecount = sectionnodecount1 + sectionnodecount2;
    
    for ielement = 1:nbeamnode-1  % combine ends # 1 and 2 of elements
        elementnodecount(ielement) = sectionnodecount(ielement) + sectionnodecount(ielement+1);
        elementsectioncount(ielement) = sectionnodecount(ielement) + sectionnodecount(ielement+1);
        elementchord(ielement) =      (sectionmaxchord(ielement) + sectionmaxchord(ielement+1))/2;
        elementchordrot(ielement) =  (sectionchordrot(ielement) + sectionchordrot(ielement+1))/2;
    end
    
    maxcountoutline = max(elementnodecount);           % determine max nodes associated with any element
    
    elementoutline = zeros(maxcountoutline,(nbeamnode-1)*3);
    
    for ielement = 1:nbeamnode-1                     % assemble single array of element outlines and sections
        matrix1 = zeros(sectionnodecount1(ielement),3);
        matrix2 = zeros(sectionnodecount2(ielement),3);
        matrix3 = zeros(sectionnodecount1(ielement+1),3);
        matrix4 = zeros(sectionnodecount2(ielement+1),3);
        
        elementnodecoord = zeros(maxcountoutline,3);         % put relevant coordinates into one array
        
        R = eye(3,3);
        if convertflag == 1
            R(1,1) = cos(yrotate(i));         % transformation matrix.
            R(1,2) = 0.0;
            R(1,3) = -sin(yrotate(i));
            R(2,1) = sin(xrotate(i))*sin(yrotate(i));
            R(2,2) = cos(xrotate(i));
            R(2,3) = sin(xrotate(i))*cos(yrotate(i));
            R(3,1) = cos(xrotate(i))*sin(yrotate(i));
            R(3,2) =-sin(xrotate(i));
            R(3,3) = cos(xrotate(i))*cos(yrotate(i));
        end
        
        for j = 1:sectionnodecount1(ielement)
            temp1(1:3) = sectionnodecoord1(ielement,j,1:3);
            temp = [xlocate(ielement);ylocate(ielement);0];
            temp2 = R*(temp1-temp);
            matrix1(j,1:3) = temp2(1:3)';
        end
        
        for j = 1:sectionnodecount2(ielement)
            temp1(1:3) = sectionnodecoord2(ielement,j,1:3);
            temp = [xlocate(ielement);ylocate(ielement);0];
            temp2 = R*(temp1-temp);
            matrix2(j,1:3) = temp2(1:3)';
        end
        
        for j = 1:sectionnodecount1(ielement+1)
            temp1(1:3) = sectionnodecoord1(ielement+1,j,1:3);
            temp = [xlocate(ielement+1);ylocate(ielement+1);0];
            temp2 = R*(temp1-temp);
            matrix3(j,1:3) = temp2(1:3)';
        end
        
        for j = 1:sectionnodecount2(ielement+1)
            temp1(1:3) = sectionnodecoord2(ielement+1,j,1:3);
            temp = [xlocate(ielement+1);ylocate(ielement+1);0];
            temp2 = R*(temp1-temp);
            matrix4(j,1:3) = temp2(1:3)';
        end
        
        %     'size of section matrices'
        [nr1,nc1] = size(matrix1);
        [nr2,nc2] = size(matrix2);
        [nr3,nc3] = size(matrix3);
        [nr4,nc4] = size(matrix4);
        elementoutline(1:nr1,(ielement-1)*3+1:(ielement-1)*3+3) =  matrix1;
        elementoutline(nr1+1:nr1+nr2,(ielement-1)*3+1:(ielement-1)*3+3) =  matrix2;
        elementoutline(nr1+nr2+1:nr1+nr2+nr3,(ielement-1)*3+1:(ielement-1)*3+3) =  matrix3;
        elementoutline(nr1+nr2+nr3+1:nr1+nr2+nr3+nr4,(ielement-1)*3+1:(ielement-1)*3+3) =  matrix4;
        
    end            % end of loop on elements
    
    [nroutline,ncoutline] = size(elementoutline);
    
end              % end of check for fea model

n1 = nbeamnode-1;
fprintf(fid4,'%5i   !_#_of_elements \n',n1);          %  print no of beam elements

%  write basic element data
for ielement = 1:nbeamnode-1;             % loop on beam elements
    
    length = sqrt((zbeamnode(ielement+1)-zbeamnode(ielement))^2 + ((xlocate(ielement+1)-xlocate(ielement))^2.0)...
        + (ylocate(ielement+1)-ylocate(ielement))^2);
    massperlen = nodemass(ielement)/length;
    
    if feaflag == 1;
        cofpxoffset = (sectioncofpx(ielement)+sectioncofpx(ielement+1))/2;    %  average of center of  pressure at all nodes
        cofpyoffset = (sectioncofpy(ielement)+sectioncofpy2(ielement+1))/2;     %  average of center of  pressure at all nodes
        LExoffset = (sectionLE(ielement,1)+sectionLE(ielement+1,1))/2;
        LEyoffset = (sectionLE(ielement,2)+sectionLE(ielement+1,2))/2;
        Ixxperlen = Ixx(ielement);
        Iyyperlen = Iyy(ielement);
        Ixyperlen = Ixy(ielement);
        cgxoffset = cgx(ielement);
        cgyoffset = cgy(ielement);
        maxchord = (sectionmaxchord(ielement)+sectionmaxchord(ielement+1))/2;
        chordrot = (sectionchordrot(ielement)+sectionchordrot(ielement+1))/2;
        Izzperlen = massmtz2(ielement)/length;
    else
        cofpyoffset = sectioncofpy(ielement);
        cofpxoffset = sectioncofpx(ielement);
        Ixxperlen = Ixx(ielement);
        Iyyperlen = Iyy(ielement);
        Izzperlen = massmtz2(ielement)/length;
        LExoffset = sectionLE(ielement,1);
        LEyoffset = sectionLE(ielement,2);
        cgxoffset = cgx(ielement);
        cgyoffset = cgy(ielement);
        maxchord = sectionmaxchord(ielement);
        chordrot = sectionchordrot(ielement);
        
    end
    
    fprintf(fid4,'%5i %9.3f %9.3f %9.2f %9.3f %9.3f %9.3f %9.3f %9.3f !elmntNo, zcoord, chord, chordRot, cgx, cgy, LEx, LEy, mass/len \n',...
        ielement,zbeamnode(ielement), maxchord,chordrot*57.3,cgxoffset,cgyoffset,...
        LExoffset,LEyoffset,massperlen);
    
    fprintf(fid5,'%5i %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n',ielement,massperlen,cofpyoffset,Ixxperlen,Iyyperlen,Izzperlen,Ixyperlen);
    fprintf(fid5,'%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n',cgxoffset,cgyoffset,maxchord,chordrot,LExoffset,LEyoffset);
    
    matrix = storesectstiff((ielement-1)*8+1:(ielement-1)*8+6,1:6);
    fprintf(fid5,'%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n',matrix);
    fprintf(fid5,'\n');
    
    
end                        % end loop on elements

fprintf(fid4,'%5i %9.3f  !end_node_no,end_node_coord \n',nbeamnode,zbeamnode(nbeamnode));          %end coordinate
fprintf(fid4,'\n');

[nrstiff,ncstiff] = size(storesectstiff);

for ielement = 1:nbeamnode-1              %  loop for element stiffness matrices etc
    
    matrix = storesectstiff((ielement-1)*8+1:(ielement-1)*8+6,1:6);
    fprintf(fid4,'%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n',matrix);
    fprintf(fid4,'\n');
    
    %                                         % calc some inertias per unit length
    
end                                              %  end of loop of elements stiffness matrices and given k & inertias

fprintf(fid4,'\n');

if feaflag ==1;                % print outline if fea model present
    for irow = 1:maxcountoutline                               %  loop on rows in section array
        fprintf(fid4,'%12.4f', elementoutline(irow,:));
        fprintf(fid4,'\n');
    end
    
end                                % end of fea model check

status=fclose(fid4);
status=fclose(fid5);

return
