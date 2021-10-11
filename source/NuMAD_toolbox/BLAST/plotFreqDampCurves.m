function plotFreqDampCurves(OmegaArray,convergedDamp,convergedFreq,analysisType)
                        
[~,numModes] = size(convergedFreq);
%plot freq vs. operating condition and damp vs. operating
%condition curves
plotMarkers={'.k','.k','.b','.b','.g','.g','.c','.c','.m','.m','ok','ok','ob','ob','og','og','oc','oc','om','om'};
figure(1);
subplot(2,1,2);
hold on;
for i=1:numModes
   plot(OmegaArray,convergedDamp(:,i),plotMarkers{i});
end
    if(strcmpi(analysisType,'F'));
        xlabel('Rotor Speed (RPM)');
    end
    if(strcmpi(analysisType,'P'));
        xlabel('Freestream Velocity');
    end
ylabel('Damping Ratio');
grid on;

subplot(2,1,1);
hold on;
for i=1:numModes
   plot(OmegaArray,convergedFreq(:,i),plotMarkers{i});
end
    if(strcmpi(analysisType,'F'));
        xlabel('Rotor Speed (RPM)');
    end
    if(strcmpi(analysisType,'P'));
        xlabel('Freestream Velocity');
    end
ylabel('Frequency (Hz)');
grid on;
end
