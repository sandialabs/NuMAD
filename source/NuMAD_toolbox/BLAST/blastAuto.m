function [convergedFreq,convergedDamp] = blastAuto(infile,outfile,OmegaArray,tol,analysisType)

%Function provides a wrapper to allow for a command line interface to a
%"struct" interface code.

%form struct from command line input option
inStruct.OmegaArray = OmegaArray;
inStruct.outFile = outfile;

a = importdata(infile);
inStruct.fstFile = a{1};
inStruct.bladeFile = a{2};
inStruct.aeroFile = a{3};

temp = load(a{4});
inStruct.pitchAxisVal = temp.pitchAxis;

blade = readFastBlade(inStruct.bladeFile);
fast = readFastMain(inStruct.fstFile);
bladeLength=fast.TurbConf.TipRad-fast.TurbConf.HubRad;
inStruct.pitchAxisDomain =  blade.prop.BlFract.*bladeLength;

inStruct.numadBladeLen = bladeLength; %no need to check for difference between
%fast and numad blade length

temp = load(a{5});
inStruct.LCS = temp.LCS;

%call to struct executive function
[convergedFreq,convergedDamp]=feaAutoAllModes(inStruct,tol,analysisType);

end

