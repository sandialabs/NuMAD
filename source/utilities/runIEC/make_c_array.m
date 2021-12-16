function make_c_array(IEC)
    figure
    numad(IEC.numadfn);
    numadinfo=guidata(gcf);
    close(gcf)
    len=length(numadinfo.station);
    carray=zeros(len,6);
    for i=(1:len)
        carray(i,1)=numadinfo.station(1,i).LocationZ; % Blade Location (m)
        carray(i,2)=(numadinfo.station(1,i).LocationZ)/(numadinfo.station(1,len).LocationZ); % BlFract
        carray(i,3)=numadinfo.station(1,i).Chord; % Chord
        carray(i,4)=numadinfo.station(1,i).Xoffset; % x-offset
        carray(i,5)=carray(i,3)-(carray(i,3)*carray(i,4));% edgewise c
        carray(i,6)=carray(i,3)*max(numadinfo.station(1,i).coords(:,2));% flapwise c
    end
    % disp('____blLoc_____BlFract___Chord_____Xoffset___edge c____flap c')
    % disp(carray)
    % figure(400)
    % plot(carray(:,2),carray(:,5),carray(:,2),carray(:,6))
    % xlabel('BlFract')
    % legend('edgewise c','flapwise c')

    save carray carray
end