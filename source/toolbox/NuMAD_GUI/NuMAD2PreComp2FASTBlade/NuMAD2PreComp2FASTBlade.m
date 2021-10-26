function varargout = NuMAD2PreComp2FASTBlade(varargin)
% NUMAD2PRECOMP2FASTBlade Convert NuMAD *.nmd and MatDBsi.txt files into PreComp
% inputs.  Then perform analysis at all blade stations and save results
% into Matlab data file PreComp_SectionData.mat and also FASTBlade.dat.
%
%  Important Notes/Limitations:
%   -?
%
%

% ========================================================
%   Written by Brian Resor, Sandia National Laboratories
%   Last update: 6/5/2013

% assume call from NuMAD
USING_GUI = true;
%    cbo = varargin{1};
data = varargin{2};
settings = varargin{3};
BATCH_RUN = varargin{4};
bmodes_path = settings.bmodes_path;
precomp_path = settings.precomp_path;
matdb = data.matdb;

if BATCH_RUN
    USING_GUI = false;
end

disp('*.nmd and MatDBsi.txt files read and blade model data loaded into workspace')

% convert information in data, Mat, and Comp to data structure for PreComp
% input creation
precomp=NuMAD2PreComp(data,matdb)

if USING_GUI, helpdlg('NuMAD-to-PreComp conversion is finished.  It is highly recommended to review all PreComp input files for consistency with NuMAD model. (Files: layup*.inp, materials.inp, shape*.inp,input*.pci)','FYI'); end
disp('NuMAD-to-PreComp conversion is finished.  It is highly recommended to review all PreComp input files for consistency with NuMAD model. (Files: layup*.inp, materials.inp, shape*.inp,input*.pci)');

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
PreComp_SectionData=runPreCompAnalysis(precomp,precomp_path,BATCH_RUN);
if BATCH_RUN
    obj = findobj(0,'Name','PreComp Analysis Results');
    close(obj);
end

% Use BModes to calculate mode shapes based on section properties from PreComp
bmodesFrequencies = BModes4NuMAD2PreComp2FASTBlade(data,bmodes_path,PreComp_SectionData.data);
if nargout > 0
    varargout{1} = bmodesFrequencies;
end

% fit polynomials to the modes
if BATCH_RUN
    modeShapes=polyfitmodes(data.batchmodelist);
else
    modeShapes=polyfitmodes;
end

% Generate FAST Blade file from this analysis
if USING_GUI || BATCH_RUN
    PreComp2FASTBlade(data,PreComp_SectionData.data,modeShapes);
    obj = findobj(0,'Name','FAST/ADAMS Blade File Parameter Inspection');
    close(obj);
else
    error('there is a problem!!')  % brr
end

disp('FAST Blade file has been written: FASTBlade_precomp.dat')

% Input file cleanup
if BATCH_RUN
    qu='Yes'; % always cleanup during batch runs
else
    qu=questdlg('Delete all miscellaneous input/output files?','File Cleanup...','Yes','No','No');
end
switch qu
    case 'Yes'
        delete('*.pci')
        delete('*.inp')
        delete PreComp_SectionData.mat
        delete polyfitdata.mat
        delete PrepMat.txt
        delete bmodes*.bmi
        delete blade_sec_props.dat
        delete bmodes*.echo
    case 'No'
end


end