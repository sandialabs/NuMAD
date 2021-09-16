function data_out = NuMAD_appdata(op,var,data_in)
%NUMAD_APPDATA  manages NuMAD application data
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   data_out = NuMAD_appdata(operation,variable,data_in)
%
%   NuMAD_appdata('init','',HMAIN) initializes the application data
%   structure. Requires HMAIN, the handle of the main gui window.
%
%   NuMAD_appdata('set',NAME,DATA) sets/creates a data field
%   called NAME in which DATA is stored
%
%   DATA = NuMAD_appdata('get',NAME) retrives the DATA stored as NAME
%
%   NuMAD_appdata('rm',NAME) removes the data field called NAME
%
%   NuMAD_appdata('cleanup') deletes all of the application data managed by
%   this function
%

instance = 1;  % for handling multiple instances of NuMAD

switch op
    case 'init'
%         if ~ishandle(data_in)
%             %jcb: during initialization, need to provide handle to main gui
%         end
        if isappdata(0, 'NuMAD_hMain')
            hMain = getappdata(0, 'NuMAD_hMain');
            hMain(instance) = data_in;  % handle to main gui of this NuMAD instance
            setappdata(0, 'NuMAD_hMain', hMain);
        else
            hMain(instance) = data_in;  % handle to main gui of this NuMAD instance
            setappdata(0, 'NuMAD_hMain', hMain);
        end
        
    case 'set'
        hMain = getappdata(0, 'NuMAD_hMain');
        setappdata(hMain(instance), var, data_in);
        
    case 'get'
        hMain = getappdata(0, 'NuMAD_hMain');
        data_out = getappdata(hMain(instance), var);
        
    case 'rm'
        hMain = getappdata(0, 'NuMAD_hMain');
        rmappdata(hMain(instance), var);
        
    case 'cleanup'
        if isappdata(0, 'NuMAD_hMain')
            hMain = getappdata(0, 'NuMAD_hMain');
            if numel(hMain)==1
                rmappdata(0, 'NuMAD_hMain');
            else
                hMain(instance)=[];
                setappdata(0, 'NuMAD_hMain', hMain);
            end
        end
    otherwise
        error('NuMAD_appdata: operation "%s" not understood',op);
end