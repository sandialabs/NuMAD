function writeCrunch(cru,output_file)
%WRITECRUNCH  Write a Crunch input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeCrunch(cru,file_name)
%    for use with Crunch v3.00.00
% 
%      cru = Crunch input structure; view default structure by utilizing 
%           readCrunch()
%      file_name = name of file to be created.

if ~exist('output_file','var')  %print to command window if no output_file given
    fid=1;
else
    fid=fopen(output_file,'wt');   %try to open output_file for Writing in Text mode
    if (fid == -1)
        error('Could not open file "%s"\n',output_file);
        return
    end
end

wip(fid,[],'-----  Crunch v3.00c Input File  -----------------------------------------------');
wip(fid,[],cru.title{1});

wip(fid,[],'-----  Job Options  ------------------------------------------------------------');
wip(fid,cru.JobOpt.Echo,'Echo:               The flag for echoing input to <root>.ech.');
wip(fid,cru.JobOpt.Out_Stats,'Out_Stats:          The flag for outputting statistics.');
wip(fid,cru.JobOpt.Out_Data,'Out_Data:           The flag for outputting modified data.');
wip(fid,cru.JobOpt.TabDelim,'TabDelim:           The flag for using tab-delimited output (best for spreadsheets).');
wip(fid,cru.JobOpt.RealFmt,'RealFmt:            The numerical-output format specifier.  See manual for limitations.');
wip(fid,cru.JobOpt.Aggregate,'Aggregate:          The flag to output aggregate-analysis files.  Use false to generate separate analysis files for each input file.');
wip(fid,cru.JobOpt.AggRoot,'AggRoot: The root name of the aggregate-analysis files, if aggregates were specified above.');

wip(fid,[],'-----  Input-Data Layout  ------------------------------------------------------');
wip(fid,cru.IptDataLay.CTRow,'CTRow:              The row with the channel titles on it (zero if no titles are available of if titles are specified below).');
wip(fid,cru.IptDataLay.CURow,'CURow:              The row with the channel units on it (zero if no units are available of if units are specified below).');
wip(fid,cru.IptDataLay.FDRow,'FDRow               The first row of data.');
wip(fid,cru.IptDataLay.NumRecs,'NumRecs:            The number of data records to read from each file (0 to automatically determine which rows to read).');
wipcsv(fid,cru.IptDataLay.TStartTEnd,'TStart, TEnd:       The start and end times (enter zeros if you want to use all the data records in the file).');

wip(fid,[],'-----  Channel Information  ----------------------------------------------------');
wip(fid,cru.ChanInfo.NumInCols,'NumInCols:          The total number of channels in each input file.');
wip(fid,cru.ChanInfo.NumCols,'NumCols:            The number channels to be used from the input files.');
wip(fid,[],'ChanTitle     ChanUnits    OrigChan #  Scale  Offset     NumCols rows of data follow.  Title and units strings must be 10 characters or less.');

Table = {cru.ChanInfo.ChanTitle,cru.ChanInfo.ChanUnits,cru.ChanInfo.OrigChan,cru.ChanInfo.Scale,cru.ChanInfo.Offset};
Frmt = {'%s    ', '%s     ', '%i    ','%6.3f    ','%6.3f    '};
wiptbl(fid,Frmt,Table,cru.ChanInfo.NumCols);

wip(fid,[],'-----  Filtering  --------------------------------------------------------------');
wip(fid,[],'0              NumFilt:            The output channels are to be modified by the IIR filter.');
wip(fid,[],'0              FiltCols:           The list of channels to filter.  Ignored if NumFilt is zero.');
wip(fid,[],'0              FiltType:           The type of filter (1-LowPass, 2-HighPass, 3-BandPass)');
wip(fid,[],'0.0            LoCut:              The low cutoff frequency (ignored for low-pass filters)');
wip(fid,[],'0.0            HiCut:              The high cutoff frequency (ignored for high-pass filters)');

wip(fid,[],'-----  Calculated Channels  ----------------------------------------------------');
wip(fid,cru.CalcChan.NumCChan,'NumCChan:           The number calculated channels to generate.');
wip(fid,cru.CalcChan.Seed,'Seed:               The integer seed for the random number generator (-2,147,483,648 to 2,147,483,647)');
wip(fid,[],'Col_Title  Units  Equation         Put each field in double quotes.  Titles and units are limited to 10 characters.  NumCChan rows of data follow.');

Table = {cru.CalcChan.ColTitle,cru.CalcChan.Units,cru.CalcChan.Equation};
Frmt = {'%s    ', '%s     ', '%s    '};
wiptbl(fid,Frmt,Table,cru.CalcChan.NumCChan);

wip(fid,[],'-----  Moving Averages  --------------------------------------------------------');
wip(fid,[],'0              NumMA:              Number of channels that will have moving averages generated for them.');
wip(fid,[],'Title   Channel #   Averaging Period');

wip(fid,[],'-----  Time and Wind Speed  ----------------------------------------------------');
wip(fid,[],'1              TimeCol:            The channel containing time.');
wip(fid,[],'0              WS_Col:             The primary wind-speed channels (used for mean wind speed and turbulence intensity, 0 for none)');

wip(fid,[],'-----  Load Roses  -------------------------------------------------------------');
wip(fid,[],'0              NumRoses:           Number of load-rose channel pairs.');
wip(fid,[],'RoseTitle   0-DegreeLoad   90-DegreeLoad   #Sectors       RoseTitle in quotes and up to 8 characters.  NumRoses rows of data follow.');

wip(fid,[],'-----  Azimuth Averages  -------------------------------------------------------');
wip(fid,[],'0              NumAACols:          Number of channels to be azimuth averaged.  Next four lines ignored if 0.');
wip(fid,[],'3              AA_Cols(:):         List of channels for azimuth averages.');
wip(fid,[],'360            NumAABins:          The number of azimuth bins.');
wip(fid,[],'2              AzimCol:            The azimuth channel.');
wip(fid,[],'True		       Out_AA:             The flag for outputting azimuth averages.');

wip(fid,[],'-----  Crosstalk  --------------------------------------------------------------');
wip(fid,[],'0              NumXT:              Number of pairs of columns that will have their crosstalk removed.');
wip(fid,[],'OutCol#1   OutCol#2   XT(1,1)   XT(1,2)   XT(2,1)   XT(2,2)         NumXT rows of data follow.');

wip(fid,[],'-----  Peak Finding  -----------------------------------------------------------');
wip(fid,[],'0              NumPFCols:          Number of output columns to be modified by the peak finder.  Next line ignored if zero.');
wip(fid,[],'               PF_Cols(:):         The list of channels to be modified by the peak finder.');

wip(fid,[],'-----  Peak and Valley Listing  ------------------------------------------------');
wip(fid,[],'0              NumPLCh:            Number of channels that will have their peaks/valleys listed to a file.  Next three lines ignored if zero.');
wip(fid,[],'1		           PL_Meth:            The method of identifying peaks (1: slope change, 2: thresholds).');
wip(fid,[],'True	         WrPLtime:           The flag to include the time in the peak-list file(s)?');
wip(fid,[],'Channel   WriteTroughs?   TroughThresh   WritePeaks?   PeakThresh        NumPLCh rows of data follow.');

wip(fid,[],'-----  Probablity Mass  --------------------------------------------------------');
wip(fid,[],'0              NumPMF:             The number of channels that will have PMFs generated for them.  Next two lines ignored if zero.');
wip(fid,[],'20             NumPMFBins:         The number of bins to use for the PMF.');
wip(fid,[],'Column#   Minimum   Maximum        If Min=Max=0, autocalculate them.  NumPMF rows of data follow.');

wip(fid,[],'-----  Rainflow Cycles  --------------------------------------------------------');
wip(fid,cru.NumRFCols,'NumRFCols:          The number of RF channels.  Next six lines ignored if zero.');
wip(fid,[],'1              RF_Per:             Number of seconds in the rainflow counting period.');
wip(fid,[],'False		       RF_Norm:            The flag for normalizing rainflow cycle counts by bin width.');
wip(fid,[],'False 	       RFZC_Blank:         The flag to specify the generation of spaces in tab-delimited output files where counts are zero.');
wip(fid,[],'0              NumRFRBins:         The number of rainflow range bins.  Use "0" to output the actual cycles instead of binned cycles.');
wip(fid,[],'0              NumRFMBins:         The number of rainflow means bins.  Use "1" to output ranges only.');
wip(fid,[],'Column#   Half-CycleMult   MaxRange   MinMean   MaxMean      MaxRange not use when NumRFRBins=0. MinMean, MaxMean not used when NumRFMBins=1.  NumRFCols rows of data follow.');

Table = {cru.ColNums,cru.HalfCycMult,zeros(cru.NumRFCols,1),zeros(cru.NumRFCols,1),zeros(cru.NumRFCols,1)};
Frmt = {'%i    ','%6.3f    ','%6.3f    ','%6.3f    ','%6.3f    '};
wiptbl(fid,Frmt,Table,cru.NumRFCols);

wip(fid,[],'-----  Extreme Events  ---------------------------------------------------------');
NumEEGrps=length(cru.EEGrps);
if NumEEGrps>0
    wip(fid,NumEEGrps,'NumEEGrps:          The number of groups of parameters that will have their extreme events recorded.');
    wip(fid,[],'GroupTitle(100 char max)   #ExtCols   ColList(#ExtCols long)   #InfCols(may be 0)   ColList(#InfCols long)     NumEEGrps rows of data follow.');
    for j=1:NumEEGrps
        wip(fid,cru.EEGrps{j},'');
    end
else
    wip(fid,[],'0              NumEEGrps:          The number of groups of parameters that will have their extreme events recorded.');
    wip(fid,[],'GroupTitle(100 char max)   #ExtCols   ColList(#ExtCols long)   #InfCols(may be 0)   ColList(#InfCols long)     NumEEGrps rows of data follow.');
end

wip(fid,[],'-----  Summary Files  ----------------------------------------------------------');
wip(fid,[],'0              NumSFCols:          The number of channels to have statistics put in separate summary files.  Next line ignored if zero.');
wip(fid,[],'0              SF_Cols(:):         The list of summary-file channels.');

wip(fid,[],'-----  Statistical Extrapolation  ----------------------------------------------');
wip(fid,[],'0              NumESCols:          The number channels for statistical extrapolation.');
wip(fid,[],'Col#   HoursToExtrapolateTo   QuantileDesired             NumESCols rows of data follow.');

wip(fid,[],'-----  Input Files  ------------------------------------------------------------');
wip(fid,cru.NumFiles,'NumFiles:           The number of input files to read. ');
for j=1:cru.NumFiles
    wip(fid,cru.InFiles{j},[]);
end

if fid~=1, fclose(fid); end
end

%==========================================================================
%===== FUNCTION DEFINITIONS ===============================================
%==========================================================================
function wip(fid,param,descrip)
% write input file parameter
if ~any(size(param)) && ~ischar(param)
    % do nothing if param = []
    % note: used ~any(size(param)) rather than isempty(param)
    %       so that unset parameters (size = [1 0]) will still
    %       get through to the following elseif statements
elseif ischar(param)
    fprintf(fid,'%-16s ',param);  %output string
elseif isfloat(param)
    if numel(param)==1
        fprintf(fid,'%-16g ',param);  %output single number
    else
        str = sprintf(' %g',param);        %create list of numbers
        str = sprintf('"%s"',str(2:end));  %quote the list of numbers
        fprintf(fid,'%-16s ',str);          %output the quoted list
    end
end
fprintf(fid,'%s\n',descrip);
end

function wiplst(fid,param,descrip)
% write input file parameter list (one per line)
for k = 1:length(param)
    fprintf(fid,'%-16s ',param{k});  %output string
    if k==1
        fprintf(fid,'%s',descrip);
    end
    fprintf(fid,'\n');
end
end

function wiptbl(fid,frmt,table,nrows)
for r=1:nrows
    for c=1:length(table)
        col = table{c};
        if iscell(col)
            param = col{r};
        else
            param = col(r);
        end
        fprintf(fid,frmt{c},param);
    end
    fprintf(fid,'\n');
end
end

function wipcsv(fid,param,descrip)
% write input file parameter list as comma separated values
str = '';
for k = 1:length(param)
    str = strcat(str,sprintf('%d,',param(k)));
end
if length(str) > 1
    str(end) = [];
end
fprintf(fid,'%-16s %s\n',str,descrip);
end

