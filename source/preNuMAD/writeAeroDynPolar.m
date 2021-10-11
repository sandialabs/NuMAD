function writeAeroDynPolar(polar,filename,aerodynVersion,multiParam)
%writeAeroDynPolar  Write an AeroDyn airfoil performance input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeAeroDynPolar(polar,'filename') 
%      Created for AeroDyn v12.58 and v13.00. Default format is 'v13'.
%
%      polar = Airfoil performance data structure; view default structure 
%          by utilizing readAirfoilData()
%      filename = name of file to be created.
%      aerodynVersion = (optional) select between 'v12' and 'v13' formats
%      multiParam = (optional) select between 'Re' and 'Ctrl' TableID
%          parameter when using the 'v12' format

if ~exist('aerodynVersion','var') || isempty(aerodynVersion)
    aerodynVersion = 'v13.0';  % flag to switch the output file format
end

if ~exist('multiParam','var') || isempty(multiParam)
    multiParam = 'Re';  % flag to switch the v12 TableID parameter
end

switch aerodynVersion
    case {'v13','v13.0','13.0'}
        ADver = 'v13';
    case {'v12','v12.0','12.0'}
        ADver = 'v12';
        % check whether the Alpha columns are compatible and consistent
        NumAlf = [polar.Tab.NumAlf];
        if ~all(NumAlf==NumAlf(1))
            error('AeroDyn v12 format requires tables to have same Alpha');
        end
        alpha = [polar.Tab.Alpha];
        alpha_diff = alpha - repmat(alpha(:,1),1,polar.NumTabs);
        if any(alpha_diff(:) > 1e-4)
            msg = 'Alpha columns appear to differ.';
            if ~isempty(filename)
                msg = sprintf('%s Using first table''s Alpha to write %s.',msg,filename);
            else
                msg = sprintf('%s Using first table''s Alpha.',msg);
            end
            warning(msg)
        end
    otherwise
        error('aerodynVersion "%s" not recognized',aerodynVersion);
end

fw = fileWriter(filename);  % create a file writer object
fw.fmt_num = '  %-8g';  % set the format strings
fw.fmt_str = '  %-8s';
fw.fmt_lst = '%-8g ';

if isequal(ADver,'v13')
    fw.text('AeroDyn airfoil file.  Compatible with AeroDyn v13.0.');
end
fw.text(polar.title{1});
fw.text(polar.title{2});
fw.num(polar.NumTabs,'Number of airfoil tables in this file');

% check whether every table has Cm data
if all(arrayfun(@(x) isfield(x,'Cm'),polar.Tab))
    UseCM = true;
else
    UseCM = false;
end

switch ADver
    case 'v13'
        for k = 1:polar.NumTabs
            fw.num(polar.Tab(k).Re/1e6  ,'Reynolds number in millions.  For efficiency, make very large if only one table.');
            fw.num(polar.Tab(k).Ctrl    ,'Control setting');
            fw.num(polar.Tab(k).AlfaStal,'Stall angle (deg)');
            fw.num(polar.Tab(k).AOL     ,'Zero lift angle of attack (deg)');
            fw.num(polar.Tab(k).CnA     ,'Cn slope for zero lift (dimensionless)');
            fw.num(polar.Tab(k).CnS     ,'Cn at stall value for positive angle of attack');
            fw.num(polar.Tab(k).CnSL    ,'Cn at stall value for negative angle of attack');
            fw.num(polar.Tab(k).AOD     ,'Angle of attack for minimum CD (deg)');
            fw.num(polar.Tab(k).Cd0     ,'Minimum CD value');
            for j = 1:polar.Tab(k).NumAlf
                fprintf(fw.fileID,'  %7.2f',polar.Tab(k).Alpha(j));
                fprintf(fw.fileID,'  %7.5f',polar.Tab(k).Cl(j));
                fprintf(fw.fileID,'  %7.5f',polar.Tab(k).Cd(j));
                if UseCM
                fprintf(fw.fileID,'  %7.5f',polar.Tab(k).Cm(j));
                end
                fprintf(fw.fileID,'\n');
            end
            fw.text('EOT');
        end
    case 'v12'
        zerolist = zeros(1,polar.NumTabs);
        if strcmpi('Ctrl',multiParam)
            fw.list([polar.Tab.Ctrl]   ,'Control setting');
        else
            fw.list([polar.Tab.Re]/1e6 ,'Reynolds number in millions');
        end
        fw.list([polar.Tab.AlfaStal],'Stall angle (deg)');
        fw.list(zerolist            ,'No longer used, enter zero');
        fw.list(zerolist            ,'No longer used, enter zero');
        fw.list(zerolist            ,'No longer used, enter zero');
        fw.list([polar.Tab.AOL]     ,'Zero lift angle of attack (deg)');
        fw.list([polar.Tab.CnA]     ,'Cn slope for zero lift (dimensionless)');
        fw.list([polar.Tab.CnS]     ,'Cn at stall value for positive angle of attack');
        fw.list([polar.Tab.CnSL]    ,'Cn at stall value for negative angle of attack');
        fw.list([polar.Tab.AOD]     ,'Angle of attack for minimum CD (deg)');
        fw.list([polar.Tab.Cd0]     ,'Minimum CD value');
        for j = 1:polar.Tab(1).NumAlf
            fprintf(fw.fileID,'  %7.2f',polar.Tab(1).Alpha(j));
            for k = 1:polar.NumTabs
                fprintf(fw.fileID,'  %7.5f',polar.Tab(k).Cl(j));
                fprintf(fw.fileID,'  %7.5f',polar.Tab(k).Cd(j));
                if UseCM
                fprintf(fw.fileID,'  %7.5f',polar.Tab(k).Cm(j));
                end
            end
            fprintf(fw.fileID,'\n');
        end
end



end