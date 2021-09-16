function BPEGUI(nmdfn)
%GUIBPE  Plotting script for determination of BPE segment edge indices
% **********************************************************************
% *           Part of the SNL Wind Turbine Analysis Toolbox            *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% **********************************************************************
%   guibpe()
%   Detailed function description describing inputs, outputs, and other necessary information.
%
%

%===== CREDITS & CHANGELOG ================================================
% 2011.06.09  brr: initial release
% yyyy.mm.dd  initials: description

[a.station, shearweb, active, ansys, a.BladeRotation, blade] = readNuMADinput(nmdfn);

N=length(a.station);
BlLength=a.station(end).LocationZ;

% get BPE segment edges
fid2 = fopen('bpe_station_ids.txt','r'); % opens basic input file
if fid2==-1
    flag=0;
    bpesta=1; % arbitrary setting
else
    string = scandm2(fid2);
[bpesta nbpesta] = sscanf(string,'%f,'); % puts values into vector
fclose(fid2);
flag=1;
end

j=1;
for i=1:N
    chord(i)=a.station(i).Chord;
    z(i)=a.station(i).LocationZ;
    xoff(i)=a.station(i).Xoffset;
    if find(bpesta==i)
        segedges(j)=a.station(i).LocationZ;
        j=j+1;
    end
end

le=chord.*xoff;
te=chord.*(xoff-1);
spacing=5; %percent of span
n=1/(spacing*0.01)+1;
guides=linspace(0,BlLength,n);

figure(1000)

subplot(2,1,1)
plot(z,zeros(N),'ro',...
    guides,zeros(n),'bv',...
    z,le,'k',...
    z,te,'k')
axis equal
legend('All Numbered NuMAD Stations')
ylabel('Chord, m')
for i=1:N
    text(z(i),0.05,num2str(i))
end
title('Blue triangles are at recommended 5% span spacing.')

if flag
    subplot(2,1,2)
    plot(segedges,zeros(nbpesta),'ro',...
        guides,zeros(n),'bv',...
        z,le,'k',...
        z,te,'k')
    axis equal
    legend('NuMAD Stations Specified in "bpe\_station\_ids.txt"')
    xlabel('Span, m')
    ylabel('Chord, m')
    for i=1:nbpesta
        text(segedges(i),0.05,num2str(bpesta(i)))
    end
else
    subplot(2,1,2)
    title('"bpe\_station\_ids.txt" file was not found')
end

end
