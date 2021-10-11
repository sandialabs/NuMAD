function varargout = readWTPerfBED(wtp_input_file,bed_file)
%READWTPERFBED  Read a WTPerf BED output file.
%   bed = readWTPerfBED('wtp_input_filename','bed_filename')
%   Returns a structure containing the WTPerf BED output data.
%
%   Some output hints:
%     bed.Element = [nElement x 1]
%     bed.Azimuth = [nSect x 1]
%     bed.'data'  = [nElement x nSect x nCases]
%   where,
%        'data'   = one of the 17 data columns of output in the BED file
%        nElement = number of aerodynamic elements used in WTPerf
%        nSect    = number of Azimuthal elements used in WTPerf
%        nCases   = number of cases (parametric or specified) analyzed by WTPerf
%
%   Plotting examples for populated bed structure.
%     Plot Cl versus rotor radius for all azimuthal sectors for Case #2,
%       plot(bed.roverR,bed.Cl(:,:,2)),xlabel('r/R'),ylabel('C_l'),legend(num2str(bed.Azimuth))
%     Plot Cl versus rotor radius averaged over all azimuthal sectors for Case #2,
%        plot(bed.roverR,mean(bed.Cl(:,:,2),2)),xlabel('r/R'),ylabel('C_l'),legend('mean C_l over all sectors')
%     Plot Cl versus Azimuth for the element at the blade tip for Case #2,
%        plot(bed.Azimuth,bed.Cl(end,:,2)),xlabel('Azimuth, deg'),ylabel('C_l')
%
%   Last update: 4/17/2013 BRR
%   Compatible with WTPerf v3.10 BED output file format

wtp=readWTPerf(wtp_input_file);
nelement = wtp.NumSeg;

if (wtp.Tilt == 0 && wtp.PreCone == 0 && wtp.Yaw == 0)
    nSect=1;
else
    nSect=wtp.NumSect;
end

% open file
fid=fopen(bed_file);
if (fid == -1)
    error('Could not open input "%s"\n',bed_file);
    return
end

% get file description from header
bed.header{1}=fgetl(fid);
bed.header{2}=fgetl(fid);
for i=1:4
    line=fgetl(fid);
end

k=1;
temp=zeros(nelement,nSect,17);
while (1)
    line=fgetl(fid);
    pat = 'Blade-element data for (?<param1>[^=]*)= (?<value1>[\d\.-]*) \w*, (?<param2>[^=]*)= (?<value2>[\d\.-]*) \w*, (?<param3>[^=]*)= (?<value3>[\d\.-]*).';
    n = regexp(line,pat,'names');
    
    var=strrep(n.param1,' ','');
    val=str2num(n.value1);
    bed=setfield(bed,var,{k},val);
    
    var=strrep(n.param2,' ','');
    val=str2num(n.value2);
    bed=setfield(bed,var,{k},val);
    
    var=strrep(n.param3,' ','');
    val=str2num(n.value3);
    bed=setfield(bed,var,{k},val);
    
    for i=1:4
        line=fgetl(fid);
    end
    
    bed.Element=(1:nelement)';
    
    for i=1:nelement
        for j=1:nSect
            if nSect>1;line=fgetl(fid);end  % there's an extra line that appears only when number of sectors is greater than 1
            temp=sscanf(line,'%g',[1,17]);
            
            bed.Azimuth(j,1)=temp(2);
            bed.LocVel(i,j,k)=temp(3);
            bed.Re(i,j,k)=temp(4);
            bed.Loss(i,j,k)=temp(5);
            bed.AxialInd(i,j,k)=temp(6);
            bed.TangInd(i,j,k)=temp(7);
            bed.AirflowAngle(i,j,k)=temp(8);
            bed.AlfaD(i,j,k)=temp(9);
            bed.Cl(i,j,k)=temp(10);
            bed.Cd(i,j,k)=temp(11);
            bed.ThrustCoef(i,j,k)=temp(12);
            bed.TorqueCoef(i,j,k)=temp(13);
            bed.PowerCoef(i,j,k)=temp(14);
            bed.ThrustperLen(i,j,k)=temp(15);
            bed.TorqueperLen(i,j,k)=temp(16);
            bed.Power(i,j,k)=temp(17);
        end
        
        line=fgetl(fid);
    end
    
    if line==-1
        break;
    end
    
    k=k+1;
end

bed.roverR=wtp.seg.RElm/wtp.RotorRad;

% close file
status=fclose(fid);

varargout{1} = bed;

end