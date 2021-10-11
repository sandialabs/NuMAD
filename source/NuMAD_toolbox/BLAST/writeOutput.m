function  writeOutput(freq,damp,phase1,phase2,imagComponentSign,sectionLoc,convergedCheck,fid,analysisType,plotFlag)
%writeOutput Writes output to file and plots mode shapes
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [freqSorted,dampSorted,modesAmpSorted,modesPhaseSorted,imagCompSignSorted] = 
%    writeOutput(freq,damp,phase1,phase2,imagComponentSign,sectionLoc,
%    convergedCheck,fid,analysisType,plotFlag)
%
%   This function excepts input of the blade configuration and results
%   information, sorts modes with respect to frequency, and writes the 
%   formatted data to file and provides mode shape plots.
%
%      input:
%      freq                = unsorted frequency list
%      damp                = unsorted damping ratio list
%      phase1              = in-phase (0 degree) mode shapes
%      phase2              = out-of-phase (90 degree) mode shape
%      imagComponentSign   = eigenvalue imaginary component sign list
%      sectionLoc          = blade section location list
%      convergedCheck      = list of booleans for converged modes
%      fid                 = output file identifier
%      analysisType        = char denoting analysis type
%      plotFlag            = flag controlling plotting (1-plot, 0-no plot)
%      
%      output:
%      none

[l,w] = size(phase1(:,:,1));

numModes = length(freq);
for i=1:1:numModes
   fprintf(fid,'MODE # %i \n\n',i);
   fprintf(fid,'Frequency: %i: \n',freq(i)); 
   fprintf(fid,'Damping %i: \n\n',damp(i));
   fprintf(fid,'0 deg Mode Shape:\n');
   fprintf(fid,'U_x          U_y          U_z          theta_x     theta_y     theta_z \n');

   if(plotFlag)
   plotModeShape(freq(i),damp(i),phase1(:,:,i),...
                 phase2(:,:,i),sectionLoc,convergedCheck(i),...
                 i,analysisType);
   end
   for(j=1:l)
        fprintf(fid,'%8.6f \t',phase1(j,:,i));
        fprintf(fid,'\n');
   end
   fprintf(fid,'\n');
   
   fprintf(fid,'90 deg Mode Shape:\n');
   fprintf(fid,'U_x          U_y          U_z          theta_x     theta_y     theta_z \n');
   for(j=1:l)
        fprintf(fid,'%8.6f \t',phase2(j,:,i));
        fprintf(fid,'\n');
   end
   
   fprintf(fid,'\n\n');
   
end

     %if plotting modes bring those with converged frequencies to the front
     if(plotFlag)
        lenCC=length(convergedCheck);
        for i=1:lenCC
            index = length(convergedCheck)+1-i; 
            if(convergedCheck(index))
                figure(index);
            end
        end
     end
end

function plotModeShape(freq,damp,phase1,phase2,sectionLoc,convergedFlag,i,analysisType)
% This function plots the mode shape of the elastic axis
%
%   input:
%   freq          = frequency of mode
%   damp          = damping ratio of mode
%   phase1        = 0 degree (in phase) mode shape
%   phase2        = 90 degree (out of phase) mode shape
%   sectionLoc    = spanwise location of structural node
%   convergedFlag = flag describing if mode is converged
%   analysisType  = char denoting analysis type
%
%   output:
%   none - No explicit function output, a plot is generated.

gridA = sectionLoc./sectionLoc(length(sectionLoc));

figure(i);
hold on;


plot(gridA,phase1(:,2),'-r',gridA,phase1(:,3),'-b',...
     gridA,phase1(:,4),'-k',gridA,phase2(:,2),'--r',...
     gridA,phase2(:,3),'--b',gridA,phase2(:,4),'--k','LineWidth',2.5);

set(gcf,'OuterPosition',[0 100 700 485])
set(gca,'Position',[.05 .11 .75 .815]);
legend('Edgewise - 0^o','Flapwise - 0^o','Torsional - 0^o','Edgewise - 90^o','Flapwise - 90^o','Torsional - 90^o',-1);
xlabel('x/L');

mTextBox = uicontrol('style','text');

if(~strcmpi(analysisType,'D'))
    set(mTextBox,'String',['k = ',num2str(freq),' Hz','                    Damping = ',num2str(damp) ],'Position',[400.0 180.0 200.0 50.0],'FontSize',12);
    if(~convergedFlag)
        mTextBox2 = uicontrol('style','text');
        set(mTextBox2,'String','NOT CONVERGED','Position',[110.0 50.0 200.0 80.0],'FontSize',20);
    end
else
    set(mTextBox,'String',['Omega_div = ',num2str(freq),' RPM'],'Position',[400.0 150.0 200.0 80.0],'FontSize',12);
end
grid on;
axis([0 1 -1 1]);
% saveas(gcf,num2str(i),'bmp');
% saveas(gcf,num2str(i),'fig');
end