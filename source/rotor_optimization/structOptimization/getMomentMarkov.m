function markov=getMomentMarkov(rccdata,wt,Yr,simtime,markovSize,chSpan,direction)
    
    if chSpan==1
        baseStr='Root';
    else
        baseStr=['Spn' int2str(chSpan-1)];
    end

    %Search through first windspeed column of rccdata for the first appearance of
    %baseStr + M + direction e.g. RootMyb1
    ct=1; %At the end of the while loop ct's value will be the row index value
          % of baseStr + M + direction e.g. RootMyb1.
    while ~(contains(rccdata{ct,1}.label,baseStr) && contains(rccdata{ct,1}.label,'M') && contains(rccdata{ct,1}.label,direction))
        ct=ct+1;
    end
    
    means=[]; %Mean fatigue data, appended for each wind speed
    ampl=[]; %Amplitude fatige data, appended for each wind speed
    cycles=[];
    for w=1:size(rccdata,2)
        % put data from rccdata structure for this channel and wind speed into a temporary variable, data
        data=rccdata{ct,w};

        % Make sure that fatigue data are only summed accross
        % windspeeds for the same channel.
        if  w==1 || strcmp(data.label,rccdata{ct,w-1}.label)
             means =[means; data.means]; %Append a column (each column of data corresponds to the data for wind speed w
             ampl =[ampl; data.amplitudes];
             cycles =[cycles; data.cycles/simtime * (60*60*24*365.24*Yr) * wt(w);];  % cycles per second times number of seconds in Yr years times Rayleigh weight factor
        else     
            error('Data channel name from current wind speed does not match the previous wind speed. Fatigue cycles cannot be summed.')
        end    
    end

    [Ni,EDGESi,BINi] = histcounts(ampl,markovSize);    %Ni, Nj, number of elements in each bin
    [Nj,EDGESj,BINj] = histcounts(means,markovSize);   %EDGES, bin edges,
                                                       %BIN, bin number assignment for each element in ampl or means

    markov=zeros(markovSize+1);
    markov(2:end,1)=0.5 * (EDGESi(1:end-1) + EDGESi(2:end))'; %Compute values at bin centers by taking the average
    markov(1,2:end)=0.5 * (EDGESj(1:end-1) + EDGESj(2:end));
    
        %indChecki=zeros(1,markovSize);
        for j=1:markovSize  %Means
            ampIndecies=find(BINj==j);
            %indCheck=0;

            for i=1:markovSize %Amplitudes
                meanIndecies=find(BINi==i);

                Indecies= intersect(ampIndecies,meanIndecies);
                markov(i+1,j+1)=sum(cycles(Indecies));
            end
        end
    %surf(markov(2:end,2:end))
end