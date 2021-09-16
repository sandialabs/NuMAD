function NuMAD2BPE2FASTBlade(varargin)
%NuMAD2BPE2FASTBlade Use ANSYS & BPE to create FAST blade file from NuMAD model.
% **********************************************************************
% *           Part of the SNL Wind Turbine Analysis Toolbox            *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% **********************************************************************
%   NuMAD2BPE2FASTBlade
%
%   Detailed function description describing inputs, outputs, and other necessary information.
%
%   nmdfn, full file name of NuMAD file describing the blade model; *.nmd
%   bmodes_path, full path to BMODES.exe.  BModes is available from NREL-NWTC.
%
%   Use BPEGUI.m to determine which NuMAD stations will be used as BPE
%   element edges in the BPE analysis.
%

%===== CREDITS & CHANGELOG ================================================
% 2011.21.09  brr: initial creation of this script
% 2011.12.22  brr: various improvements for integration with NuMAD
% yyyy.mm.dd  initials: description

if ischar(varargin{1})
    % assume filename argument
    USING_GUI = false;
    nmdfn = varargin{1};
    bmodes_path = varargin{2};
    [data.station, shearweb, active, ansys, data.BladeRotation, blade] = readNuMADinput(nmdfn);
    
    if ispc
        % get the windows %APPDATA% path,
        % ( typically C:\Users\{username}\AppData\Roaming )
        userpath = fullfile(getenv('APPDATA'),'NuMAD');
    elseif isunix
        % get the unix/linux $HOME path
        userpath = fullfile(getenv('HOME'),'.numad');
    elseif ismac
        errordlg('Mac is not currently supported by NuMAD.','Error');
        error('Mac is not currently supported by NuMAD.');
    else
        errordlg('Your system is not supported by NuMAD.','Error');
        error('Your system is not supported by NuMAD.');
    end
    settings = readNuMADsettings(fullfile(userpath,'settings.txt'));
else
    % assume call from NuMAD
    USING_GUI = true;
    cbo = varargin{1};
    data = varargin{2};
    settings = varargin{3};
    bmodes_path = settings.bmodes_path;
    
    callapp = guidata(cbo);
    app.caller = callapp.fh;  % store figure handle provided by calling script
    app.numadpath = callapp.numadpath;
    app.userpath = callapp.userpath;
end

if isempty(bmodes_path)
    if USING_GUI, errordlg('Path to BModes not specified. Aborting.','Error: BModes Path'); end
    error('Path to BModes not specified. Aborting.');
end

disp('*.nmd file read and blade model data loaded into workspace')

flag=exist('bpe_station_ids.txt','file');
if ~flag
    if USING_GUI, warndlg('Please define the BPE segments first. ("bpe_station_ids.txt" not found)','Notification'); end
    return;
end

flag=exist('displacement.txt','file');
if ~flag  % Generate files that feed into displacement.txt by performing ANSYS static analyses
    
    %jcb: should we check that the boundary condition is cantilevered?
    writeAPDL4BPE
    
    if isequal(0,exist('shell7bpe.src','file'))
        if USING_GUI, warndlg('Generate an ANSYS "shell7.src" file before beginning BPE analysis.','Notification'); end
        return;
    end
    
    % run ANSYS
    msgg='BPE requires nodal displacements to be computed by ANSYS.  Computed displacements are stored in the file, "displacement.txt."  This file is not found from a previous ANSYS analysis.  Executing new ANSYS batch run using shell7bpe.src in order to generate a new "displacement.txt" file for subsequent BPE analysis.';
    disp(msgg)
    if USING_GUI, helpdlg(msgg,'Notification'); end
    
    if isempty(settings.ansys_path)
        if USING_GUI, errordlg('Path to ANSYS not specified. Aborting.','Error: ANSYS Call'); end
        error('Path to ANSYS not specified. Aborting.');
    end
    
    [status,result] = dos(sprintf('"%s" -b -p %s -I shell7bpe.src -o output.txt',settings.ansys_path,settings.ansys_product));
    
    if isequal(status,0)
        % dos command completed successfully; log written to output.txt
        if USING_GUI, helpdlg('ANSYS batch run has completed. See "output.txt" for any warnings.','ANSYS Call Completed'); end
        disp('ANSYS batch run has completed. See "output.txt" for any warnings.');
    end
    
    if ~isequal(status,0)
        % an error has occured which is stored in output.txt
        if USING_GUI, errordlg('Could not complete ANSYS call. See "output.txt" for details.','Error: ANSYS Call'); end
        error('Could not complete ANSYS call. See "output.txt" for details.');
    end
    
    % combine individual saved results into a single file, displacement.txt
    disp('Combining ANSYS results data files into single file: displacement.txt')
    one=dlmread('nodeloc.txt');
    two=dlmread('fdisp.txt');
    thr=dlmread('mdisp.txt');
    two=two/1e3;
    thr=thr/1e3;
    fou=dlmread('nMasses.txt');
    fiv=[one two thr fou];
    file4='displacement.txt';
    fid=fopen(file4,'wt');
    if (fid == -1)
        errordlg(sprintf('Could not open file "%s"',file4),'Error');
        return
    end
    fprintf(fid,'%i\t%10.8e\t%10.8e\t%10.8e\t%10.8e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%18.16e\t%10.8e\n',fiv');
    fclose(fid);
    
    delete('nodeloc.txt');
    delete('fdisp.txt');
    delete('mdisp.txt');
    delete('nMasses.txt');
else
    disp('The "displacement.txt" file is found -- Skipping ANSYS analysis of the blade.');
end

% get BPE segment edges
fid2 = fopen('bpe_station_ids.txt','r'); % opens basic input file
if fid2==-1, error('BPE segment edges input file "bpe_station_ids.txt" not found'), end
string = scandm2(fid2);
[bpesta nbpesta] = sscanf(string,'%f,'); % puts values into vector

j=1;
segedges=[];
tw_aero=[];
for i=1:length(data.station)
    if find(bpesta==i)
        segedges(j)=data.station(i).LocationZ;
        tw_aero(j)=data.station(i).DegreesTwist;
        j=j+1;
    end
end
precurve=zeros(1,length(segedges));
save numaddata data segedges tw_aero

% write information to bpe input file, bpe.txt
writeBPEInput(segedges,precurve);
disp('BPE input file written: bpe.txt')

ebeam131a % requires bpe.txt and displacement.txt
disp('BPE analysis complete')

BPEPost % requires numadguidata and bpe.txt
disp('Postprocessing of BPE data complete.  Data saved in BPE_SectionData.mat')

BModes4NuMAD2BPE2FASTBlade(bmodes_path)
disp('BMODES analysis of beam properties complete.  Data saved in bmodes.bmi, bmodes_sec_prop.dat & bmodes.out')

polyfitmodes % requires bmodes.bmi and bmodes.out
disp('Polynomial fit of modeshapes complete')

if exist('nmdfn','var')
    BPEPost2FASTBlade(nmdfn) % requires BPE_SectionData and polyfitdata
else
    BPEPost2FASTBlade(data.station) % requires BPE_SectionData and polyfitdata
end

disp('FAST Blade file has been written: FASTBlade_bpe.dat')

qu=questdlg('Delete all miscellaneous input/output files?','File Cleanup...','Yes','No','No');

switch qu
    case 'Yes'
        delete bmodes.bmi
        delete bmodes_sec_props.dat
        delete bmodes.out
        delete bmodes.echo
        delete *.asv
        delete bpe.txt
        delete kout.txt
        delete xlsout.txt
        delete admout.txt
        delete eigstore.txt
        delete flexmatrixsym.txt
        delete sectflexmatrix.txt
        delete sectiondispl.txt
        delete sectstiffmatrix.txt
%         delete singularout.txt
        delete stiffmatrixsym.txt
%         delete displacement.txt
        delete polyfitdata.mat
        delete numaddata.mat
    case 'No'
end

end