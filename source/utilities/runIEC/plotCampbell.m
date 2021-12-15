function plotCampbell(IEC)
    global fastPath
    global adamsPath  %not used?
    global turbsimPath
    global iecwindPath
    global crunchPath
    global mbcPath
    % initialize arrays
    rotCampbell=[];
    rotSpd=[];
    tmp=[];
    id=1;l=0;

    % point to MBC analysis scripts
    addpath(mbcPath);

    % if ~isfield(IEC,'ratedSpeed')
    if 0
        IEC.ratedSpeed=1;
        IEC.nondimFlag=false;
    else
        IEC.nondimFlag=true;
    end

    %% Analyze system modes

    % for ws=[0 IEC.lin]
    for ws=[0 (5:2:22)]
        save vars ws IEC rotCampbell rotSpd id l tmp
        clear all
        load vars

        fn=sprintf('out\\Linear_%dmps.lin',ws);
        RootName=fn(1:end-4);
        GetMats
        if 0
            mbc3
            inFreqs=MBC_NaturalFrequencyHz;
            inMag=MBC_ModeShapeMagnitude;
            inPha=MBC_ModeShapePhaseDeg;
        else
            cce
            inFreqs=AvgNaturalFrequencyHz;
            inMag=AvgModeShapeMagnitude;
            inPha=AvgModeShapePhaseDeg;
        end

        [freq,ix]=sort(inFreqs);  % sort frequencies
        mag=inMag(:,ix);  % sort based on same order as sorted frequencies
        pha=inPha(:,ix);  % sort based on same order as sorted frequencies

        tmp{id}.freq=freq;
        l=max([l length(freq)]);
        rotSpd=[rotSpd RotSpeed*60/2/pi];

        if ws==0 && 0   % plot degrees of freedom for parked rotor
            n=(1:size(mag,1));
            for i=1:size(inMag,2)

                figure(2)
                subplot(4,5,i)
                plot(n,mag(:,i),'b-o')
                title(sprintf('%4.1f RPM - #%i - %4.2f Hz',RotSpeed*60/2/pi,i,tmp{id}.freq(i)))
                xlabel('DOF')
                ylabel('Mag')
                grid on
                set(gcf,'Position',[260 130 1400 850]);

                figure(3)
                subplot(4,5,i)
                plot(n,pha(:,i),'kv')
                title(sprintf('%4.1f RPM - #%i - %4.2f Hz',RotSpeed*60/2/pi,i,tmp{id}.freq(i)))
                xlabel('DOF')
                ylabel('Phase')
                grid on
                set(gcf,'Position',[260 130 1400 850]);

            end

    %         return;
        end

        id=id+1;

    end

    DescStates
    rmpath(mbcPath);

    for i=1:length(tmp)
        temp=[zeros(l-length(tmp{i}.freq),1); tmp{i}.freq];
        rotCampbell=[rotCampbell temp];
    end

    %% Plot system modes

    f=60/IEC.ratedSpeed;  % rated rotor rotating frequency in 1/Hz
    f=1;

    figure(4)
    plot(rotSpd/IEC.ratedSpeed,rotCampbell*f,'-x',...
        rotSpd/IEC.ratedSpeed,rotSpd*1/60*f,'k-',...
        rotSpd/IEC.ratedSpeed,rotSpd*3/60*f,'k-',...
        rotSpd/IEC.ratedSpeed,rotSpd*6/60*f,'k-',...
        rotSpd/IEC.ratedSpeed,rotSpd*1.1/60*f,'k:',...
        rotSpd/IEC.ratedSpeed,rotSpd*3.3/60*f,'k:',...
        rotSpd/IEC.ratedSpeed,rotSpd*6.3/60*f,'k:',...
        rotSpd/IEC.ratedSpeed,rotSpd*0.9/60*f,'k:',...
        rotSpd/IEC.ratedSpeed,rotSpd*2.7/60*f,'k:',...
        rotSpd/IEC.ratedSpeed,rotSpd*5.7/60*f,'k:')
    grid on

    if IEC.nondimFlag
    xlabel({'\Omega_R/\Omega_{R,rated}' sprintf('(\\Omega_{R,rated}=%4.1f rpm)',IEC.ratedSpeed)})
    %     ylabel('f/\Omega_{R,rated}')
        ylabel('Natural frequency (Hz)')
    else
        xlabel('\Omega_R (rpm)')
        ylabel('f (Hz)')
    end

    title({'Full system aeroelastic frequencies','Computed and analyzed by FAST/AeroDyn and MBC'})
    set(gcf,'Position',[680 343 560 755]);
    axis([0 5 0 15])
    YTick = (0:1:25);
    YTickLabel = YTick;
    set(gca,'YTick',YTick);
    set(gca,'YTickLabel',YTickLabel);


    %% Plot Blade modes

    bldCampbell=[];
    bldSpd=(0:3:15);

    for omega=bldSpd
        fn=sprintf('NuMAD\\bmodes_%drpm.out',omega);
        bmout=readBModesOut(fn,13);
        bldCampbell=[bldCampbell bmout.freq'];
    end

    figure(5)
    plot(bldSpd/IEC.ratedSpeed,bldCampbell'*f,'-x',...
        bldSpd/IEC.ratedSpeed,bldSpd*1/60*f,'k-',...
        bldSpd/IEC.ratedSpeed,bldSpd*3/60*f,'k-',...
        bldSpd/IEC.ratedSpeed,bldSpd*6/60*f,'k-',...
        bldSpd/IEC.ratedSpeed,bldSpd*1.1/60*f,'k:',...
        bldSpd/IEC.ratedSpeed,bldSpd*3.3/60*f,'k:',...
        bldSpd/IEC.ratedSpeed,bldSpd*6.3/60*f,'k:',...
        bldSpd/IEC.ratedSpeed,bldSpd*0.9/60*f,'k:',...
        bldSpd/IEC.ratedSpeed,bldSpd*2.7/60*f,'k:',...
        bldSpd/IEC.ratedSpeed,bldSpd*5.7/60*f,'k:')
    grid on

    if IEC.nondimFlag
        xlabel(sprintf('Non-dimensional rotor speed (RPM_r_a_t_e_d=%4.1f)',IEC.ratedSpeed))
        ylabel('Non-dimensional natural frequency')
    else
        xlabel('Rotor speed (RPM)')
        ylabel('Natural frequency (Hz)')
    end

    title({'Campbell diagram','Single blade aeroelastic frequencies','Computed and analyzed by BModes'})

    save Campbell bldCampbell rotCampbell rotSpd bldSpd
end