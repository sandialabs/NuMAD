function [output] = readDLCoutput(filename)
%readDLCoutput   Read in the runIEC output file, 'IECDLC_Results.slsx' to 
%                convert it into a structure within the MATLAB workspace.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
[num, txt, raw] = xlsread(filename);
% find the variable headers in the data file
varIndex = [1; find(isnan(num(:,1))==1)+1]; 
% number of rows corresponding to the number of Design Load Cases tested
numRows = varIndex(2)-varIndex(1)-1;

for ii = 1:length(varIndex)
    % set up the string name as a structure cell containing a table
    strVar = txt{varIndex(ii),2};
    output.(strVar) = table();
    % save the columns of the data table within the structure cell
    output.(strVar).name = txt(varIndex(ii)+1:varIndex(ii)+numRows,1);
    output.(strVar).data = num(varIndex(ii):varIndex(ii)+numRows-1,1);
    output.(strVar).time = num(varIndex(ii):varIndex(ii)+numRows-1,2);
    output.(strVar).chan = txt(varIndex(ii)+1:varIndex(ii)+numRows,4);
    output.(strVar).file = txt(varIndex(ii)+1:varIndex(ii)+numRows,5);
end

