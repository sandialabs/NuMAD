function NuMAD2PreComp2FASTBlade(varargin)
% NUMAD2PRECOMP2FASTBlade Convert NuMAD *.nmd and MatDBsi.txt files into PreComp
% inputs.  Then perform analysis at all blade stations and save results
% into Matlab data file PreComp_SectionData.mat and also FASTBlade.dat.
%

if ischar(varargin{1})
    % assume filename argument
    USING_GUI = false;
    nmdfn = varargin{1};
    bmodes_path = varargin{2};
    precomp_path = varargin{3};

    [data.station, data.shearweb, data.active, ansys, data.BladeRotation, blade] = readNuMADinput(nmdfn);
    matdb = readMatDB('MatDBsi.txt');
else
    % assume call from NuMAD
    USING_GUI = true;
%    cbo = varargin{1};
    data = varargin{2};
    settings = varargin{3};
    bmodes_path = settings.bmodes_path;
    precomp_path = settings.precomp_path;
    matdb = data.matdb;
end

disp('*.nmd and MatDBsi.txt files read and blade model data loaded into workspace')

% convert information in data, Mat, and Comp to data structure for PreComp
% input creation
precomp=NuMAD2PreComp(data,matdb);
helpdlg('NuMAD-to-PreComp conversion is finished.  It is highly recommended to review all PreComp input files for consistency with NuMAD model. (Files: layup*.inp, materials.inp, shape*.inp,input*.pci)','FYI')

if isempty(precomp_path)
    if USING_GUI, errordlg('Path to PreComp not specified. Aborting.','Error: PreComp Path'); end
    error('Path to PreComp not specified. Aborting.');
end

if isempty(bmodes_path)
    if USING_GUI, errordlg('Path to BModes not specified. Aborting.','Error: BModes Path'); end
    error('Path to BModes not specified. Aborting.');
end

% run PreComp analyses at each of the stations and then combine results
% from each analysis into one Matlab array.  Save as PreComp_SectionData.mat
runPreCompAnalysis(precomp, precomp_path);

% Use BModes to calculate mode shapes based on section properties from PreComp
BModes4NuMAD2PreComp2FASTBlade(data,bmodes_path);

% fit polynomials to the modes
polyfitmodes

% Generate FAST Blade file from this analysis
if USING_GUI
    PreComp2FASTBlade(data.station);
else
    PreComp2FASTBlade(nmdfn)
end

disp('FAST Blade file has been written: FASTBlade_precomp.dat')

% Input file cleanup
qu=questdlg('Delete all miscellaneous input/output files?','File Cleanup...','Yes','No','No');
switch qu
    case 'Yes'
        delete('*.pci')
        delete('*.inp')
        delete PreComp_SectionData.mat
        delete polyfitdata.mat
        delete bmodes.out
        delete PrepMat.txt
        delete bmodes.bmi
        delete blade_sec_props.dat
        delete bmodes.echo
    case 'No'
end


end