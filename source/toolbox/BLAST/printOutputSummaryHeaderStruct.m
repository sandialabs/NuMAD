function [fid] = printOutputSummaryHeaderStruct(outfile,inStruct,analysisType)
%printOutputSummaryHeader  prints header for BLAST output summary file
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [fid] = printOutputSummaryHeader(outfile,inputfile,analysisType)
%                    
%   This function prints header information for the BLAST output summary
%   file including analysis date/time, input file name, and analysis type.
%
%      input:
%      outfile       = root filename for output files
%      inputfile     = input filename for analysis
%      analysisType  = character denoting analysis type
 
%      output:
%      fid           = file id for output summary (.sum) file

fid = fopen([outfile,'.sum'],'wt');
fprintf(fid,'Sandia National Laboratories BLade Aeroelastic Stability Tool (BLAST)\n');
clckinfo = clock;
fprintf(fid,'%s \t %2.2d:%2.2d\n\n',date,clckinfo(4),clckinfo(5));

fprintf(fid,'Working Directory:\t\t%s\n\n',pwd);
fprintf(fid,'FAST File:\t\t%s\n',inStruct.fstFile);
fprintf(fid,'Blade File:\t\t%s\n',inStruct.bladeFile);
fprintf(fid,'AeroDyn File:\t\t%s\n\n',inStruct.aeroFile);

if(strcmp(analysisType,'F'))
    atypestr = 'Rotating flutter analysis';
end
if(strcmp(analysisType,'P'))
    atypestr = 'Parked flutter analysis';
end
fprintf(fid,'Analysis Type:\t%s\n\n',atypestr);

fprintf(fid,'Frequency and damping values written to %s\n\n',[outfile,'.mat']);

end

