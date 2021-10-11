xlsFile='SNL61p5m_j.xlsx';
%  =[panel thickness 1 2; spar cap thickness 1 2; spar cap width; te reinf];
% lb = [ 20;20;    60;60;  450;  20 ];
% ub = [ 80;80;   200;200; 700; 150 ];
lb = [  20;20;    60;  60;  450;   50 ];
ub = [  60;60;   220; 220;  650;  150 ];
A=[ -1  1  0  0  0  0;
    0  0  0  0  0  0;
    0  0  0  0  0  0;
    0  0  0  0  0  0;
    0  0  0  0  0  0;
    0  0  0  0  0  0 ];
b=zeros(length(A),1);
pop=10;
gen=50;
restartFile='layupDesignCandidates.txt';
% restartFile='';
layupDesignOpt(xlsFile,lb,ub,A,b,pop,gen,restartFile)

% layupDesignCheck(xlsFile)

N=12; numelem=2.^linspace(1,8,N); aeSizes=numelem.^(-0.5);
% aeSizes=[0.05 0.1 0.15 0.20 0.25 0.3 0.4]';
% layupMeshCheck(xlsFile,aeSizes,0)

% layupDesignIterate(xlsFile)
