function writePreCompInput(precomp,output_file)
%WRITEPRECOMPINPUT  Write a PreComp main input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writePreCompInput(precomp,'output_file')
%      Created for PreComp v1.00.03
%      Designed for use with NuMAD2PreComp()
%
%      precomp = PreComp data structure; view default structure 
%          in NuMAD2PreComp()
%      output_file = filename of file to be created.
%

for i=1:length(precomp.station)
    dot=max(strfind(output_file,'.'));
    ext=output_file(dot+1:end);
    fn=output_file(1:dot-1);
    num_output_file=[fn '_' num2str(i) '.' ext];
    
    if ~exist('num_output_file','var')  %print to command window if no output_file given
        fid=1;
    else
        fid=fopen(num_output_file,'wt');   %try to open output_file for Writing in Text mode
        if (fid == -1)
            error('Could not open file "%s"\n',num_output_file);
            return
        end
    end
    
    
    wip(fid,[],'*****************  main input file for PreComp *****************************')
    wip(fid,[],'Sample Composite Blade Section Properties')
    wip(fid,[],'')
    wip(fid,[],'General information -----------------------------------------------')
    wip(fid,1,'Bl_length   : blade length (m)')
    wip(fid,2,'N_sections  : no of blade sections (-)')
    wip(fid,length(precomp.material.id),'N_materials : no of materials listed in the materials table (material.inp)')
    wip(fid,3,'Out_format  : output file   (1: general format, 2: BModes-format, 3: both)')
    wip(fid,'f','TabDelim     (true: tab-delimited table; false: space-delimited table)')
    wip(fid,[],'')
    wip(fid,[],'Blade-sections-specific data --------------------------------------')
    wip(fid,[],'Sec span     l.e.     chord   aerodynamic   af_shape    int str layup')
    wip(fid,[],'location   position   length    twist         file          file')
    wip(fid,[],'Span_loc    Le_loc    Chord    Tw_aero   Af_shape_file  Int_str_file')
    wip(fid,[],'  (-)        (-)       (m)    (degrees)       (-)           (-)')
    wip(fid,[],'')
    
    for j=0:1
        PropTable = {j, precomp.station(i).lepos, precomp.station(i).chlen, ...
            precomp.station(i).twist, {['shape_' num2str(i) '.inp']}, {['layup_' num2str(i) '.inp']} };
        Frmt = {'  %3.2f     ', '%4.3f   ', '%5.3e   ', ...
            '%4.2f      ', '%s     ', '%s'};
        
        wiptbl(fid,Frmt,PropTable,1);
    end
    
    wip(fid,[],'')
    wip(fid,[],'Webs (spars) data  --------------------------------------------------')
    wip(fid,[],'')
    
    try
        N_webs=length(precomp.station(i).sw.stackid);
    catch
        N_webs=0;
    end
    
    N_fb=precomp.station(i).fb.flag;
    wip(fid,N_webs+N_fb,'Nweb        : number of webs (-)  ! enter 0 if the blade has no webs')
    wip(fid,1,'Ib_sp_stn   : blade station number where inner-most end of webs is located (-)')
    wip(fid,2,'Ob_sp_stn   : blade station number where outer-most end of webs is located (-)')
    wip(fid,[],'')
    wip(fid,[],'Web_num   Inb_end_ch_loc   Oub_end_ch_loc (fraction of chord length)')
    
    for j=1:N_webs+N_fb
if j<=N_webs
    PropTable = [j precomp.station(i).sw.chloc(j) precomp.station(i).sw.chloc(j)];
else
%     PropTable = [j 1.0 1.0];
    PropTable = [j 0.9999 0.9999];
    % BRR: the line above is trying to prevent this error:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       main input file read successfully
        %  BLADE STATION   1 analysis begins ----
        %
        %  NOTE: blunt te identified between airfoil nodes  53and  54
        %  ERROR** web no 3 outside airfoil boundary at the current blade station
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
Frmt = '%5.4f        ';
        wiphlst(fid,Frmt,PropTable);
    end
    
    if fid~=1, fclose(fid); end
end
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

function wiphlst(fid,Frmt,param)
% write input file horizontal list of numbers (all on one line, separated by spaces)
for k = 1:length(param)
    fprintf(fid,Frmt,param(k));
end
fprintf(fid,'\n');
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