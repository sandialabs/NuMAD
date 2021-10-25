function make_c_array_BladeDef(blade)

len=length(blade.ispan);
carray=zeros(len,6);
for i=(1:len)
% %     carray(i,1)=numadinfo.station(1,i).LocationZ; % Blade Location (m)
% %     carray(i,2)=(numadinfo.station(1,i).LocationZ)/(numadinfo.station(1,len).LocationZ); % BlFract
% %     carray(i,3)=numadinfo.station(1,i).Chord; % Chord
% %     carray(i,4)=numadinfo.station(1,i).Xoffset; % x-offset
    carray(i,1)=blade.ispan(i); % Blade Location (m)
    carray(i,2)=blade.ispan(i)/blade.ispan(end); % BlFract
    carray(i,3)=blade.ichord(i); % Chord
    carray(i,4)=blade.xoffset(i); % x-offset
    % calculation of distance to neutral axis -- assumes a symmetric
    % structure (which is a WRONG assumption)
    carray(i,5)=carray(i,3)-(carray(i,3)*carray(i,4));% edgewise c
    carray(i,6)=0.5*blade.ipercentthick(i)/100*blade.ichord(i);% flapwise c
end
% disp('____blLoc_____BlFract___Chord_____Xoffset___edge c____flap c')
% disp(carray)
% figure(400)
% plot(carray(:,2),carray(:,5),carray(:,2),carray(:,6))
% xlabel('BlFract')
% legend('edgewise c','flapwise c')

save carray carray
end