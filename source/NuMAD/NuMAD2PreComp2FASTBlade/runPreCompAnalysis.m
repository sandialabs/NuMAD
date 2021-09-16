function runPreCompAnalysis(precomp, precomp_path)

% Perform separate analysis at each station
for i=1:length(precomp.station)
    disp(' ')
    disp(['********************** Analyzing section #' num2str(i) ' **********************'])
    dos(sprintf('"%s" input_%d.pci',precomp_path,i));
    disp(' ')
end

% Read in PreComp results and gather in one array
for i=1:length(precomp.station)
    fn=['input_' num2str(i) '.out_gen'];
    [a,ffn,nh,SR,hl,fpos]=txt2mat(fn);
    data(i,:)= [precomp.station(i).span a(1,2:end)];
    
end

% Clean up files
delete('*bmd')
delete('*gen')

% PreComp output labels and units
labels={ 'span_loc' 'chord' 'tw_aero' 'ei_flap' 'ei_lag' 'gj' 'ea' 's_fl' 's_af' 's_al' 's_ft' 's_lt' 's_at' 'x_sc' 'y_sc' 'x_tc' 'y_tc' 'mass' 'flap_iner' 'lag_iner' 'tw_iner' 'x_cm' 'y_cm'};
units={'(-)' '(m)' '(deg)' '(Nm^2)' '(Nm^2)' '(Nm^2)' '(N)' '(Nm^2)' '(Nm)' '(Nm)' '(Nm^2)' '(Nm^2)' '(Nm)' '(m)' '(m)' '(m)' '(m)' '(Kg/m)' ' (Kg-m)' ' (Kg-m)' '(deg)' '(m)' '(m)'};

% plot PreComp Analysis results
figure('Name','PreComp Analysis Results');
for i=2:size(data,2)
    subplot(5,5,i)
    plot(data(:,1),data(:,i),'b-o')
    xlabel(strrep(labels{1},'_','\_'))
    ylabel(strrep(labels{i},'_','\_'))
end

PreComp_SectionData=data;
save PreComp_SectionData PreComp_SectionData labels units

end
