function iecwind=makeIECWind(params)

global iecwindPath

fst=readFastMain([params.fstfn '.fst']);
ad=readFastAD(strrep(fst.ADFile,'"',''));

hm=pwd;
% % cd([params.parDir 'wind'])
if(~exist('wind','dir'))
    mkdir('wind');
end
cd('wind')

% populate IECWind information
iecwind.SI='True';
iecwind.TStart=params.delay+10;
iecwind.Class=params.Class;
iecwind.Turb=params.TurbClass;
iecwind.Slope=0;
iecwind.Shear=3;
iecwind.HubHt= fst.TurbConf.TowerHt+fst.TurbConf.Twr2Shft+sind(-1*fst.TurbConf.ShftTilt)*-1*fst.TurbConf.OverHang;
iecwind.RotDia=fst.TurbConf.TipRad*1.1;
iecwind.CutIn=params.operatingPoints(1);
iecwind.Rated=params.operatingPoints(2);
iecwind.CutOut=params.operatingPoints(3);
iecwind.Conditions={'ECD+R-2.0',...
    'ECD+R',...
    'ECD+R+2.0',...
    'ECD-R-2.0',...
    'ECD-R',...
    'ECD-R+2.0'};
ws=(params.operatingPoints(1):1:params.operatingPoints(3));
for i=1:length(ws)
    iecwind.Conditions{end+1}=sprintf('EWSV+%3.1f',ws(i));
    iecwind.Conditions{end+1}=sprintf('EWSV-%3.1f',ws(i));
end
iecwind.Conditions{end+1}='EWM50';
iecwind.Conditions{end+1}='EWM01';

% ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% for parallel simulation of the runIEC script ??
% ** NOTE: IECWIND WILL NOT WORK FOR INPUT FILES UNLESS NAMED IEC.IPT
ctr = 0;
inputFile = 'IEC.ipt';
while exist(inputFile,'file') ~= 0
    ctr = ctr+1;
    inputFile = ['IEC' num2str(ctr) '.ipt'];
end
writeIECWind(iecwind, 'IEC.ipt');%inputFile); % NOTE: IECWIND WILL NOT WORK FOR INPUT FILES UNLESS THEY ARE NAMES 
dos([iecwindPath ' ' inputFile],'-echo');
% delete(inputFile)
% ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

switch params.Class
    case 1
        V50=1.4*50; % m/s, average wind speed of IEC Class I site (Vref=50m/s); IEC Section 6.3.1.1 Eqn (9)
    case 2
        V50=1.4*42.5; % m/s, average wind speed of IEC Class I site (Vref=42.5m/s); IEC Section 6.3.1.1 Eqn (9)
    case 3
        V50=1.4*37.5; % m/s, average wind speed of IEC Class I site (Vref=37.5m/s); IEC Section 6.3.1.1 Eqn (9)
end
V1=0.8*V50;
shr=0.11;
time_delay = 30; % time delay to reach maximum speed [sec]

% write EWM files with yaw misalignment
for yawAngle=-180:5:175   % write 50-year files
    fn=sprintf('EWM50%+03i.wnd',yawAngle);
    disp(sprintf('Generating non-IECWind file %s',fn));
    fid=fopen(fn,'wt');
    fprintf(fid,'!----------------------------------------------------------\n');
    fprintf(fid,'! Time  Wind    Wind    Vertical  Horiz.	Pwr.Law	    Lin.Vert.   Gust\n');
    fprintf(fid,'!  	 Speed	 Dir     Speed     Shear    Vert.Shr    Shear       Speed\n');
    fprintf(fid,'! (sec) (m/s)   (deg)   (m/s)				                        (m/s)\n');
    fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',0,0,yawAngle,0,0,shr,0,0);
    fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',time_delay,V50,yawAngle,0,0,shr,0,0);
    fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',999,V50,yawAngle,0,0,shr,0,0);
    fclose(fid);
end
for yawAngle=-30:5:30   % write 1-year files
    fn=sprintf('EWM01%+03i.wnd',yawAngle);
    disp(sprintf('Generating non-IECWind file %s',fn));
    fid=fopen(fn,'wt');
    fprintf(fid,'!----------------------------------------------------------\n');
    fprintf(fid,'! Time  Wind    Wind    Vertical  Horiz.	Pwr.Law	    Lin.Vert.   Gust\n');
    fprintf(fid,'!  	 Speed	 Dir     Speed     Shear    Vert.Shr    Shear       Speed\n');
    fprintf(fid,'! (sec) (m/s)   (deg)   (m/s)				                        (m/s)\n');
    fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',0,0,yawAngle,0,0,shr,0,0);
    fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',time_delay,V1,yawAngle,0,0,shr,0,0);
    fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',999,V1,yawAngle,0,0,shr,0,0);
    fclose(fid);
end
% wind ramp file
fn='ramp.wnd';
disp(sprintf('Generating non-IECWind file %s',fn));
fid=fopen(fn,'wt');
fprintf(fid,'!----------------------------------------------------------\n');
fprintf(fid,'! Time  Wind    Wind    Vertical  Horiz.	Pwr.Law	    Lin.Vert.   Gust\n');
fprintf(fid,'!  	 Speed	 Dir     Speed     Shear    Vert.Shr    Shear       Speed\n');
fprintf(fid,'! (sec) (m/s)   (deg)   (m/s)				                        (m/s)\n');
fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',0,  params.operatingPoints(1),0,0,0,0,0,0);
fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',30, params.operatingPoints(1),0,0,0,0,0,0);
fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',999,params.operatingPoints(3),0,0,0,0,0,0);
fclose(fid);
% steady wind files
for windSpeed=[0 params.lin]
    fn=sprintf('steady_%i.wnd',windSpeed);
    disp(sprintf('Generating non-IECWind file %s',fn));
    fid=fopen(fn,'wt');
    fprintf(fid,'!----------------------------------------------------------\n');
    fprintf(fid,'! Time  Wind    Wind    Vertical  Horiz.	Pwr.Law	    Lin.Vert.   Gust\n');
    fprintf(fid,'!  	 Speed	 Dir     Speed     Shear    Vert.Shr    Shear       Speed\n');
    fprintf(fid,'! (sec) (m/s)   (deg)   (m/s)				                        (m/s)\n');
    fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',0,  windSpeed,0,0,0,0,0,0);
    fprintf(fid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n',999,windSpeed,0,0,0,0,0,0);
    fclose(fid);
end

cd(hm)

iecwind = {};
end