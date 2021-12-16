function xformedcoords = xform(coords,offset,chord,twist)
% xformedcoords = xform(coords,offset,chord,twist)
%   Transform normalized airfoil into actual blade geometry.
% Created 1/26/2011 BRR

SS1=coords;
SS1(:,1)=SS1(:,1)-offset;
SS1=SS1*chord;
t=twist;
rot=[cos(t) -sin(t);
    sin(t) cos(t)];

for j=1:length(SS1)
    SS1(j,:)=SS1(j,:)*rot;
end

xformedcoords=SS1;

end
