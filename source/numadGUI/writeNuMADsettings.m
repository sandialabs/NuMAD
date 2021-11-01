function writeNuMADsettings(settings,filename)
%WRITENUMADSETTINGS  Output the NuMAD settings file
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeNuMADsettings(settings,FILENAME)
%     Save NuMAD settings out to file (typically 'settings.txt'
%     located in %APPDATA%\NuMAD\)
%
%   See also readNuMADsettings

%===== CREDITS & CHANGELOG ================================================
%Developed by Wind & Water Power Technologies, Sandia National Laboratories
%2011.02.16  JCB: first draft based on writeNuMADinput()

% ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% don't write settings file when in parallel
workID = getCurrentWorker;
if isempty(workID) % not running on parallel node/worker - write settings
    
    % Open the file for Writing in Text mode
    fid = fopen(filename,'wt');
    if (fid == -1)
        error('Could not open file "%s"',filename);
    end 
    
    try
        ReleaseVersion = 'v3.0.0';
        fprintf(fid,'<numad_settings release="%s">\n',ReleaseVersion);
        fprintf(fid,'  <ansys_path>%s</ansys_path>\n',settings.ansys_path);
        fprintf(fid,'  <ansys_product>%s</ansys_product>\n',settings.ansys_product);
        fprintf(fid,'  <bmodes_path>%s</bmodes_path>\n',settings.bmodes_path);
        fprintf(fid,'  <precomp_path>%s</precomp_path>\n',settings.precomp_path);
        fprintf(fid,'  <gui_layout monitors="%s">\n',settings.monitors);
        fprintf(fid,'    <xy_main>%d %d %d %d</xy_main>\n',settings.xy_main);
        fprintf(fid,'    <xy_materials>%d %d</xy_materials>\n',settings.xy_materials);
        fprintf(fid,'    <xy_isotropic>%d %d</xy_isotropic>\n',settings.xy_isotropic);
        fprintf(fid,'    <xy_orthotropic>%d %d</xy_orthotropic>\n',settings.xy_orthotropic);
        fprintf(fid,'    <xy_composite>%d %d</xy_composite>\n',settings.xy_composite);
        fprintf(fid,'    <xy_activelist>%d %d</xy_activelist>\n',settings.xy_activelist);
        fprintf(fid,'    <xy_ansys>%d %d</xy_ansys>\n',settings.xy_ansys);
        fprintf(fid,'    <xy_genline>%d %d</xy_genline>\n',settings.xy_genline);
        fprintf(fid,'    <xy_plot3d>%d %d</xy_plot3d>\n',settings.xy_plot3d);
        fprintf(fid,'    <xy_flutterinput>%d %d</xy_flutterinput>\n',settings.xy_flutterinput);
        fprintf(fid,'    <xy_pathconfig>%d %d</xy_pathconfig>\n',settings.xy_pathconfig);
        fprintf(fid,'    <xy_bpesegments>%d %d %d %d</xy_bpesegments>\n',settings.xy_bpesegments);
        fprintf(fid,'  <gui_layout>\n');
        fprintf(fid,'  <job_path>%s</job_path>\n',settings.job_path);
        fprintf(fid,'</numad_settings>\n');
        
        % Close the file
        fclose(fid);
    catch ME
        % Close the file
        fclose(fid);
        % rethrow error
        rethrow(ME);
    end
    
else % running on parallel worker -- don't write settings file
    % do nothing.
    disp('Parallel Operation: do not write NuMAD settings file...')
end
% ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
