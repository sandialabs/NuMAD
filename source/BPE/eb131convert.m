function [xyData] = eb131convert(nodeData,sectionData,storeData)
%
%  converts the stiffness matrices from global to "local" coordinates.
%  13 april 2004
%  22 july 2005, ref to profile removed
%

nbeamnode = nodeData.nbeamnode;
zbeamnode = nodeData.zbeamnode;
sectionmaxchord = sectionData.sectionmaxchord;
sectionchordrot = sectionData.sectionchordrot;
sectionLE = sectionData.sectionLE;
storestiffsym = storeData.storestiffsym;
storesectstiff = storeData.storesectstiff;

%
xlocate = zeros(nbeamnode,1);
ylocate = zeros(nbeamnode,1);
xrotate = zeros(nbeamnode,1);
yrotate = zeros(nbeamnode,1);
matrix = zeros(6,6);
localstiffsym = zeros((nbeamnode-1)*8,6);
localsectstiff = zeros((nbeamnode-1)*8,6);
E = zeros(6,6);
E(1,5) = 1;
E(2,4) = -1;

coordoffset = 0.50;       % set "center" at 0.5 of chord
%                          %  loop on sections
for isection = 1:nbeamnode
    
    %  determine x,y coords of curved section based on fraction of chord
    xlocate(isection) = sectionLE(isection,1) + sectionmaxchord(isection)*coordoffset*sin(sectionchordrot(isection));
    ylocate(isection) = sectionLE(isection,2) + sectionmaxchord(isection)*coordoffset*cos(sectionchordrot(isection));
end                    %  end of loop on section
%
for ielement = 1:nbeamnode-1                        %  loop on beam elements
    xoff = (xlocate(ielement+1)+xlocate(ielement))/2;    % calc mean offset to middle of element
    yoff = (ylocate(ielement+1)+ylocate(ielement))/2;
    i1 = (ielement-1)*8;
    matrix  = storesectstiff(i1+1:i1+6,1:6);   % extract section stiffness wrt global coords
    
    for i=1:6                                         % partially transform matrix for translation
        matrix(i,4) = matrix(i,4)-yoff*matrix(i,3);
        matrix(i,5) = matrix(i,5)+xoff*matrix(i,3);
        matrix(i,6) = matrix(i,6)+yoff*matrix(i,1)-xoff*matrix(i,2);
    end
    
    for j=1:6                                         % finish transforming matrix for translation
        matrix(4,j) = matrix(4,j)-yoff*matrix(3,j);
        matrix(5,j) = matrix(5,j)+xoff*matrix(3,j);
        matrix(6,j) = matrix(6,j)+yoff*matrix(1,j)-xoff*matrix(2,j);
    end
    %  calc rotation of curved center line
    xrotate(ielement) = atan((ylocate(ielement+1) - ylocate(ielement))/(zbeamnode(ielement+1)-zbeamnode(ielement)));
    yrotate(ielement) =-atan((xlocate(ielement+1) - xlocate(ielement))/(zbeamnode(ielement+1)-zbeamnode(ielement)));
    
    R = zeros(3,3);
    R(1,1) = cos(yrotate(isection));         % transformation matrix
    R(1,2) = 0.0;
    R(1,3) = -sin(yrotate(isection));
    R(2,1) = sin(xrotate(isection))*sin(yrotate(isection));
    R(2,2) = cos(xrotate(isection));
    R(2,3) = sin(xrotate(isection))*cos(yrotate(isection));
    R(3,1) = cos(xrotate(isection))*sin(yrotate(isection));
    R(3,2) =-sin(xrotate(isection));
    R(3,3) = cos(xrotate(isection))*cos(yrotate(isection));
    
    Rzero = zeros(3,3);
    R66 = [R,Rzero;Rzero,R];
    
    matrix  = R66*matrix*R66';              % transform to local coords
    storesectstiff(i1+1:i1+6,1:6) = matrix;
    
    length = ((zbeamnode(ielement+1)-zbeamnode(ielement))^2 + ((xlocate(ielement+1)-xlocate(ielement))^2)...
        + (ylocate(ielement+1)-ylocate(ielement))^2 )^0.5;
    vector = [1 1 1 1 1 1]*length;           % create H matrix
    H = diag(vector);
    H(4,2) = -(length^2)/2;
    H(5,1) = - H(4,2);
    Q = diag(vector)*length/2;               % create Q matrix
    Q(4,2) = -(length^3)/3;
    Q(5,1) = -Q(4,2);
    %                                           % form element stiffness matrix
    matrix = inv(inv(matrix)*H + E*inv(matrix)*Q);
    storestiffsym(i1+1:i1+6,1:6) = matrix;
    %
    
end                                 % end of loop on elements

xyData.xlocate = xlocate;
xyData.ylocate = ylocate;
xyData.xrotate = xrotate;
xyData.yrotate = yrotate;

end





