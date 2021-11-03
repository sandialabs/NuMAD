function PreComp_SectionData = runPreCompAnalysis(precomp,precomp_path,batch_run)

% Perform separate analysis at each station
for i=1:length(precomp.station)
    if 0 % display outputs
    disp(['********************** Analyzing section #' num2str(i) ' **********************'])
    [s,w] = dos(sprintf('"%s" input_%d.pci',precomp_path,i))
    disp(' ')
    else
        if 1 % display results of loading file to screen
            disp(['********************** Analyzing section #' num2str(i) ' **********************'])
            disp(sprintf('"%s" input_%d.pci',precomp_path,i))
            [s,w] = dos(sprintf('"%s" input_%d.pci',precomp_path,i))
        else
            str='dos(sprintf(''"%s" input_%d.pci'',precomp_path,i))';
%             evalc(str); % ble: changed to a direct dos command because I
%             was having problems with this section freezing.
            dos(sprintf('"%s" input_%d.pci',precomp_path,i))
        end
    end
%     pause
end

% Read in PreComp results and gather in one array
for i=1:length(precomp.station)
    fn=['input_' num2str(i) '.out_gen'];
    [a,ffn,nh,SR,hl,fpos]=txt2mat(fn,'InfoLevel',0);
    try%ble
    data(i,:)= [precomp.station(i).span a(1,2:end)];
    catch%ble
        keyboard
    end
end

% Clean up files
delete('*bmd')
delete('*gen')

% PreComp output labels and units
labels={ 'span_loc' 'chord' 'tw_aero' 'ei_flap' 'ei_lag' 'gj' 'ea' 's_fl' ...
    's_af' 's_al' 's_ft' 's_lt' 's_at' 'x_sc' 'y_sc' 'x_tc' 'y_tc' 'mass' ...
    'flap_iner' 'lag_iner' 'tw_iner' 'x_cm' 'y_cm'};
units={'(-)' '(m)' '(deg)' '(Nm^2)' '(Nm^2)' '(Nm^2)' '(N)' '(Nm^2)' '(Nm)' ...
    '(Nm)' '(Nm^2)' '(Nm^2)' '(Nm)' '(m)' '(m)' '(m)' '(m)' '(Kg/m)' ...
    ' (Kg-m)' ' (Kg-m)' '(deg)' '(m)' '(m)'};

% plot PreComp Analysis results
if ~batch_run
    figure('Name','PreComp Analysis Results');
    for i=2:size(data,2)
        subplot(5,5,i)
        plot(data(:,1),data(:,i),'b-o')
        xlabel(strrep(labels{1},'_','\_'))
        ylabel(strrep(labels{i},'_','\_'))
    end
end

PreComp_SectionData.data = data;
PreComp_SectionData.labels = labels;
PreComp_SectionData.units = units;

end
