%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Part of the SNL NuMAD Toolbox                    
%  Developed by Sandia National Laboratories Wind Energy Technologies 
%              See license.txt for disclaimer information             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef PolarDef < handle 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ``PolarDef``  A class definition for airfoil polar data.
%
% Examples: 
% 
%	``pol = PolarDef();``
% 
%	``pol = PolarDef(FILENAME);``
%
% Where FILENAME is the file containing airfoil polar data
%
%	``pol.rawlist``  is the data column names
% 
%	``pol.param``    is groups of Key=Value parameter structures
%
% See also ``AirfoilDef``, ``xlsBlade``, ``BladeDef``, ``StationDef``
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        file = struct('name','','date','','bytes',0);  % stucture describing data file
        source                 % source of data (cell array of strings)
        titleLine              % title line from file
        notes                  % notes from file; user may add more cells manually
        param = struct();      % grouped parameters (source, AeroDyn, ...)
        rawlist                % data column names for rawdata
        rawdata                % polar data, raw from source
        modopts = cell(0);     % options used to modify polar data
        modlist                % data column names for moddata
        moddata                % polar data, modified
    end
    
    properties (Hidden)
        fileLinesCellArray
    end
    
    methods
        % Constructor for PolarDef objects
        function obj = PolarDef(filedescrip,fileformat)
            
            % Initialize the "intrisic" properties structure
            obj.param.intrinsic = struct(...
                'maxCL',0,...
                'maxCLAlpha',0,...
                'zeroCLAlpha',0,...
                'maxLoD',0,...
                'maxLoDLift',0,...
                'maxLoDAlpha',0,...
                'alphaTrend',[0 10],...
                'SlopeIntercept',[0 0],...
                'manual',false);
            if nargin > 0
                if isstruct(filedescrip)
                    % filedescrip can be stucture, or ...
                    if isfield(filedescrip,'name')
                        obj.file.name = filedescrip.name;
                    else
                        error('PolarDef:fileDescripError','input Structure not recognized');
                    end
                    if isfield(filedescrip,'date')
                        obj.file.date = filedescrip.date;
                    end
                    if isfield(filedescrip,'bytes')
                        obj.file.bytes = filedescrip.bytes;
                    end
                elseif ischar(filedescrip)
                    % ... filedescrip can be filename string
                    obj.file.name = filedescrip;
                end
                
                % Open the file and read the entire contents
                fcontents = freadcontents(obj.file.name);
                obj.fileLinesCellArray = regexp(fcontents,'[^\n]+','match');
                
                % Read the other file attributes if filedescrip was string
                if ischar(filedescrip)
                    tmp = dir(filedescrip);
                    obj.file.date = tmp.date;
                    obj.file.bytes = tmp.bytes;
                end
                
                % Determine which read method to use
                if (nargin > 1) && ~isempty(fileformat)
                    switch upper(fileformat)
                        case {'XFOIL','XFOIL 6.96','XFOIL 6.97'}
                            obj.readXFoil(obj.fileLinesCellArray);
                        case {'CACTUS'}
                            obj.readCactus(obj.fileLinesCellArray);
                        case {'AERODYN'}
                            obj.readAerodyn(obj.fileLinesCellArray);
                        otherwise
                            error('PolarDef:formatNotRecognized','Specified file format not recognized: %s',fileformat);
                    end
                elseif ~isempty(strfind(fcontents,'XFOIL'))
                    if ~isempty(strfind(fcontents,'Version 6.96')) || ...
                       ~isempty(strfind(fcontents,'Version 6.97'))
                        obj.readXFoil(obj.fileLinesCellArray);
                    else
                        warning('Support for this XFOIL version not tested.\n  polar file "%s"\n',obj.file.name);
                        obj.readXFoil(obj.fileLinesCellArray);
                    end
                else
                    error('PolarDef:fileNotRecognized','File format not recognized: %s',obj.file.name);
                end
                                
                % Prepare initial modlist (could be followed by 
                % a call to obj.resetModData)
                obj.modlist = {'alpha','CL','CD','CM'};
                obj.resetModData;
                obj.updateIntrinsic;
            end
        end
                
        function x = getRawData(obj,xvar)
            % ``x = obj.getRawData(VARNAME);``
            % where ``'VARNAME'`` is a string from obj.rawlist
            % and matching is not case sensitive ('CL'=='cl')
            %
            % Example:
            %
            %	``x = obj.getRawData('alpha');``
            %
            msg = sprintf(['Method getRawData requires scalar object:\n',...
                           '   obj(k).getRawData   if numel(obj) > 1']);
            assert(isscalar(obj),msg)
            
            indDiv = strfind(xvar,'/');
            if isempty(indDiv)
                xind = strcmpi(xvar,obj.rawlist);  % not case sensitive
                assert(any(xind),'variable "%s" not recognized',xvar);
                x = obj.rawdata(:,xind);
            else
                % recursive call until there is no '/'
                % using indDiv(end) preserves left-to-right order of 
                %   operation, while indDiv(1) does not
                xvar1 = xvar(1:(indDiv(end)-1)    );
                xvar2 = xvar(  (indDiv(end)+1):end);
                x = obj.getRawData(xvar1) ./ obj.getRawData(xvar2);
            end
            
        end
        
        function x = getModData(obj,xvar)
            % ``x = obj.getModData(VARNAME);``
            % where ``'VARNAME'`` is a string from obj.modlist
            % and matching is not case sensitive ('CL'=='cl')
            %
            % Example:
            %
            %	``x = obj.getModData('alpha');``
            %
            msg = sprintf(['Method getModData requires scalar object:\n',...
                           '   obj(k).getModData   if numel(obj) > 1']);
            assert(isscalar(obj),msg)
            
            indDiv = strfind(xvar,'/');
            if isempty(indDiv)
                xind = strcmpi(xvar,obj.modlist);  % not case sensitive
                assert(any(xind),'modlist label "%s" not recognized',xvar);
                x = obj.moddata(:,xind);
            else
                % recursive call until there is no '/'
                % using indDiv(end) preserves left-to-right order of 
                %   operation, while indDiv(1) does not
                xvar1 = xvar(1:(indDiv(end)-1)    );
                xvar2 = xvar(  (indDiv(end)+1):end);
                x = obj.getModData(xvar1) ./ obj.getModData(xvar2);
            end
        end
        
        function x = getParam(obj,group,label)
            Nobj = numel(obj);
            x = zeros(Nobj,1);
            for k=1:Nobj
                if isfield(obj(k).param,group)
                    if isfield(obj(k).param.(group),label)
                        pval = obj(k).param.(group).(label);
                    else
                        error('Parameter label "%s" not found',label);
                    end
                else
                    error('Parameter group "%s" not found',group);
                end
                assert(isscalar(pval),'Only scalar parameters allowed.');
                x(k) = pval;
            end
        end
        
        function xind = matchModList(obj,cell_xvar)
            msg = sprintf(['Method matchModList requires scalar object:\n',...
                           '   obj(k).matchModList   if numel(obj) > 1']);
            assert(isscalar(obj),msg)
            N = numel(cell_xvar);
            xind = zeros(1,N);
            for k=1:N
                tmp = find(1==strcmpi(cell_xvar{k},obj.modlist));  % not case sensitive
                assert(~isempty(tmp),'modlist label "%s" not recognized',cell_xvar{k});
                assert(numel(tmp)==1,'modlist label "%s" found multiple times',cell_xvar{k});
                xind(k) = tmp;
            end
        end
        
        function plotRaw(obj,xvar,yvar,varargin)
            % ``obj.plotRaw(XVAR,YVAR);``
            % where ``'XVAR'`` and ``'YVAR'`` are strings from
            % ``obj.rawlist``
            %
            % Example:
            %
            %	``obj.plotRaw('alpha','CL');``
            %            
            N = numel(obj);
            for k = 1:N
                x = obj(k).getRawData(xvar);
                y = obj(k).getRawData(yvar);
                if k>1, hold on; end  % "add" subsequent plots
                plot(x,y,varargin{:})
                if isequal(k,1),holdState=ishold; end; % remember original hold state
            end
            if ~holdState
                hold off; % reset original hold state if "off"
            end
%             obj.recolorplot(1:N,get(gca,'ColorOrder'));
        end
        
        function plotMod(obj,xvar,yvar,varargin)
            % ``obj.plotMod(XVAR,YVAR);``
            % where ``'XVAR'`` and ``'YVAR'`` are strings from
            % ``obj.modlist``
            %
            % Example:
            %
            %	``obj.plotMod('alpha','CL');``
            %
            N = numel(obj);
            for k = 1:N
                x = obj(k).getModData(xvar);
                y = obj(k).getModData(yvar);
                if k>1, hold on; end  % "add" subsequent plots
                plot(x,y,varargin{:})
                if isequal(k,1),holdState=ishold; end; % remember original hold state
            end
            if ~holdState
                hold off; % reset original hold state if "off"
            end
%             obj.recolorplot(1:N,get(gca,'ColorOrder'));
        end
        
        function resetModData(obj)
            Nobj = numel(obj);
            for j=1:Nobj
                Ncols = numel(obj(j).modlist);
                Nrows = size(obj(j).rawdata,1);
                table = nan(Nrows,Ncols);
                for k=1:Ncols
                    varstr = obj(j).modlist{k};
                    try
                        table(:,k) = obj(j).getRawData(varstr);
                    catch
                        warning('resetModData skipping modlist{''%s''}\n',varstr);
                    end
                end      
                obj(j).moddata = table;
            end
        end
        
        function updateModData(obj)
            obj.resetModData;
            Nobj = numel(obj);
            for j=1:Nobj
                Nmod = numel(obj(j).modopts);
                for k=1:Nmod
                    mopts = obj(j).getModOpts(k);
                    switch upper(mopts.ModType)
                        case 'UPDATEINTRINSIC'
                            obj(j).updateIntrinsic;
                        case '3DSTALL'
                            obj(j).apply3DStall(mopts);
                        case 'EXTRAP'
                            obj(j).applyExtrap(mopts);
                        case 'DYNSTALL'
                            obj(j).calcDynStall(mopts);
                        case 'RESAMPLE'
                            obj(j).applyResample(mopts);
                    end
                end
            end
        end
        
        function mopts = getModOpts(obj,k)
            msg = sprintf(['Method getModOpt requires scalar object:\n',...
                           '   obj(j).getModOpt(k)   if numel(obj) > 1']);
            assert(isscalar(obj),msg)
            mopts = obj.modopts{k};
            switch upper(mopts.ModType)
                case '3DSTALL'
                    if ~isfield(mopts,'AlphaTrend')
                        mopts.AlphaTrend = obj.param.intrinsic.alphaTrend;
                    end
                case 'EXTRAP'
                    
                case 'DYNSTALL'
                    if ~isfield(mopts,'AlphaTrend')
                        mopts.AlphaTrend = obj.param.intrinsic.alphaTrend;
                    end
                    if ~isfield(mopts,'StallAngle')
                        mopts.StallAngle = obj.param.intrinsic.maxCLAlpha;
                    end
                case 'RESAMPLE'
                    
            end
        end
        
        function addModOpts(obj,varargin)
            % OPTIONS for ``'3DStall'``:
            %
            %	``ModType`` = ``'3DStall'``
            %
            %	``RPM`` = rotor speed in rpm
            %
            %	``R`` = rotor radius R
            %
            %	``V`` = wind speed
            %
            %	``rOverR`` = r/R location of airfoil section
            %
            %	``Chord`` = chord at r/R
            %
            %	``AlphaEnd`` = end of correction
            %
            %	``AlphaTrend`` = [minAlpha maxAlpha] for CL slope calc
            %
            % OPTIONS for ``'Extrap'``:
            %
            %	``ModType`` = ``'Extrap'``
            %
            %	``CDMax`` = maximum drag coefficient in extrapolation
            %
            %	``UseCM`` = use CM data (if available)
            %
            %	``AlphaFP`` = (optional) Get CD from flat-plate theory for abs(angle) > AlphaFP. Viterna method (recommended) used when AlphaFP not set or empty [].
            %
            % OPTIONS for ``'DynStall'``:
            %
            %	``ModType`` = ``'DynStall'``
            %
            %	``StallAngle`` = stall angle in degrees
            %
            %	``NegStallCn`` = Cn at stall for negative angle of attack
            %
            %   ``AlphaTrend`` = ``[minAlpha maxAlpha]`` for CN slope calc
            %
            % OPTIONS for ``'Resample'``:
            %
            %	``ModType`` = ``'Resample'``
            %
            %	``Alpha`` = array of Alpha values
            %
            %	``method`` = interpolation method:
            %	``'linear'``, ``'spline'``, ``'pchip'``
            %
            mopts = PolarDef.createModOpts(varargin{:});
            N = numel(obj);
            for k=1:N
                obj(k).modopts{end+1} = mopts;
            end
        end
        
        function clearModOpts(obj)
            N = numel(obj);
            for k=1:N
                obj(k).modopts = cell(0);
            end
        end
        
        function apply3DStall(obj,mopts)
            % ``obj.apply3DStall(mopts)``   
            %
            % Example:
            %
            %	``mopts =
            %	struct('RPM',12,'R',35,'V',8,'rOverR',0.9,'Chord',1.2,'AlphaEnd',14,'AlphaTrend',[0 5])``
            %
            RPM      = mopts.RPM;
            R        = mopts.R;
            V        = mopts.V;
            rOverR   = mopts.rOverR;
            Chord    = mopts.Chord;
            AlphaEnd = mopts.AlphaEnd;
            AlphaMinTrend = min(mopts.AlphaTrend);
            AlphaMaxTrend = max(mopts.AlphaTrend);
            if isfield(mopts,'abd')
                % Basic constants of Selig method.  Do not change except for research
                % default a=1, b=1, d=1   (abd = [1 1 1])
                abd = mopts.abd;
            else
                % use defaults
                abd = [];
            end
            
            N = numel(obj);
            for k=1:N
                columns = obj(k).matchModList({'alpha','CL','CD'});
                Table2D = obj(k).moddata(:,columns);
                Table3D = Selig_3DStall(Table2D,RPM,R,V,rOverR,Chord,AlphaMinTrend,AlphaMaxTrend,AlphaEnd,abd);
                [~,ia] = setxor(Table2D(:,1),Table3D(:,1));
                obj(k).moddata(ia,:) = [];  % remove rows missing from Table3D
                obj(k).moddata(:,columns) = Table3D;
            end
        end
        
        function applyExtrap(obj,mopts)
            CDMax  = mopts.CDMax;
            UseCM  = mopts.UseCM;
            
            if isfield(mopts,'AlphaFP')
                % get CD from flat plate theory for abs(angles) > AlphaFP
                AlphaFP = mopts.AlphaFP;
            else
                % use Viterna method for CD
                AlphaFP = [];
            end
            
            if UseCM
                matchList = {'alpha','CL','CD','CM'};
            else
                matchList = {'alpha','CL','CD'};
            end
            
            N = numel(obj);
            for k=1:N
                columns = obj(k).matchModList(matchList);
                InputTable = obj(k).moddata(:,columns);
                OutputTable = TableExtrap(InputTable,CDMax,UseCM,AlphaFP);
                [nrows,ncols] = size(OutputTable);
                obj(k).moddata = nan(nrows,ncols);
                obj(k).moddata(:,columns) = OutputTable;
                %jcb: need to handle case when size(moddata,2) > numel(columns)
            end
        end
        
        function applyResample(obj,mopts)
            Alpha = mopts.Alpha;
            method = mopts.method;
            
            N = numel(obj);
            for k=1:N
                column = obj(k).matchModList({'alpha'});
                InputTable = obj(k).moddata;
                OutputTable = interp1(InputTable(:,column),InputTable,Alpha(:),method);
                obj(k).moddata = OutputTable;
            end
        end
        
        function calcDynStall(obj,mopts)
            StallAngle = mopts.StallAngle;
            NegStallCn = mopts.NegStallCn;
            AlphaMinTrend = min(mopts.AlphaTrend);
            AlphaMaxTrend = max(mopts.AlphaTrend);
            
            matchList = {'alpha','CL','CD'};
            N = numel(obj);
            for k=1:N
                columns = obj(k).matchModList(matchList);
                InputTable = obj(k).moddata(:,columns);
                obj(k).param.AeroDyn = DynStall(InputTable,AlphaMinTrend,AlphaMaxTrend,StallAngle,NegStallCn);
            end
        end
        
        function updateIntrinsic(obj)
            Nobj = numel(obj);
            for j=1:Nobj
                if isequal(true,obj(j).param.intrinsic.manual)
                    continue;
                end
                if isempty(obj(j).moddata)
                    alpha = obj(j).getRawData('alpha');
                    CL = obj(j).getRawData('CL');
                    CD = obj(j).getRawData('CD');
                else
                    alpha = obj(j).getModData('alpha');
                    CL = obj(j).getModData('CL');
                    CD = obj(j).getModData('CD');
                end
                [maxCL,maxCLAlpha] = obj(j).findMaxCL(alpha,CL);
                obj(j).param.intrinsic.maxCL = maxCL;
                obj(j).param.intrinsic.maxCLAlpha = maxCLAlpha;
                obj(j).param.intrinsic.zeroCLAlpha = obj(j).findZeroCL(alpha,CL);
                [alphaTrend,SlopeIntercept] = obj(j).findAlphaTrend(alpha,CL);
                obj(j).param.intrinsic.alphaTrend = alphaTrend;
                obj(j).param.intrinsic.SlopeIntercept = SlopeIntercept;
                [maxLoD,k] = max(CL./CD);
                obj(j).param.intrinsic.maxLoD = maxLoD;
                obj(j).param.intrinsic.maxLoDLift = CL(k);
                obj(j).param.intrinsic.maxLoDAlpha = alpha(k);
            end
        end
        
        function setIntrinsic(obj,varargin)
            if nargin<=1
                obj.param.intrinsic.manual = false;
                return
            end
            for k=2:2:(nargin-1)
                label = varargin{k-1};
                value = varargin{k};
                if ischar(label) && isfield(obj.param.intrinsic,label)
                    if isscalar(value) && isfloat(value)
                        obj.param.intrinsic.(label) = value;
                        obj.param.intrinsic.manual = true;
                    end
                else
                    error('unknown field "%s"',label);
                end
            end
        end
        
        function afC = interp1(obj,X,Xq,alpha1,alpha3)
            
            % X:  interpolation variable (one for each airfoil)
            % Xq: query points at which polars are to be interpolated
            if nargin<=3
                alpha1 = [];
            end
            if nargin<=4
                alpha3 = [];
            end
            Nobj = numel(obj);
            [X,Xind] = unique(X);
            assert(length(X)==Nobj,'X must have unique value for each airfoil');
            assert(all(Xq>=X(1)),'desired Xq out of bounds of X');
            assert(all(Xq<=X(Nobj)),'desired Xq out of bounds of X');
            Nxq = numel(Xq);
            afC(Nxq) = PolarDef; % pre-allocate array
            for kxq = 1:Nxq
                xquery = Xq(kxq);
                kb = find(X>xquery,1,'first');
                if isempty(kb)
                    kb = length(X);
                end
                afA = obj(Xind(kb-1));
                afB = obj(Xind(kb  ));
                w = (xquery - X(kb-1))/(X(kb  ) - X(kb-1));
                afC(kxq) = afA.blend(afB,w,alpha1,alpha3);
            end
        end
        
        function af = blend(self,other,w,alpha1,alpha3)
            if nargin<=3 || isempty(alpha1)
                alpha1 = 0;
            end
            if nargin<=4 || isempty(alpha3)
                alpha3 = 12;
            end
            if isempty(self.moddata)
                alphaA = self.getRawData('alpha');
                CLA = self.getRawData('CL');
                CDA = self.getRawData('CD');
                CMA = self.getRawData('CM');
            else
                alphaA = self.getModData('alpha');
                CLA = self.getModData('CL');
                CDA = self.getModData('CD');
                CMA = self.getModData('CM');
            end
            if isempty(other.moddata)
                alphaB = other.getRawData('alpha');
                CLB = other.getRawData('CL');
                CDB = other.getRawData('CD');
                CMB = other.getRawData('CM');
            else
                alphaB = other.getModData('alpha');
                CLB = other.getModData('CL');
                CDB = other.getModData('CD');
                CMB = other.getModData('CM');
            end
            minAlpha = max(alphaA(1),alphaB(1));
            maxAlpha = min(alphaA(end),alphaB(end));
            
            NA = numel(alphaA);
            indA = zeros(NA,1);
%             alpha2 = self.param.intrinsic.maxCLAlpha;
            alpha2 = self.param.intrinsic.maxLoDAlpha;
            for k = 1:NA
                if alphaA(k) < alpha1
                    % create index which runs from -2 to -1 between
                    % minAlpha and alpha1
                    indA(k) = -1 + (alphaA(k)-alpha1)/(alpha1-minAlpha);
                    if alphaA(k) <= minAlpha
                        kA1a = k;
                    end
                    kA1b = k;
                elseif alphaA(k) < alpha2
                    % create index which runs from -1 to 0 between
                    % alpha1 and alpha2
                    indA(k) = (alphaA(k) - alpha2)/(alpha2-alpha1);
                    kA2 = k;
                elseif alphaA(k) <= alpha3
                    % create index which runs from 0 to 1 between
                    % alpha2 and alpha3
                    indA(k) = (alphaA(k) - alpha2)/(alpha3-alpha2);
                    kA3a = k+1;
                else
                    % create index which runs from 1 to 2 between
                    % alpha3 and maxAlpha
                    indA(k) = 1 + (alphaA(k)-alpha3)/(maxAlpha-alpha3);
                    if alphaA(k) <= maxAlpha
                        kA3b = k;
                    end
                end
            end
            if kA3a >= NA, kA3a=NA; kA3b=NA-1; end;
            
            NB = numel(alphaB);
            indB = zeros(NB,1);
%             alpha2 = other.param.intrinsic.maxCLAlpha;
            alpha2 = other.param.intrinsic.maxLoDAlpha;
            for k = 1:NB
                if alphaB(k) < alpha1
                    % create index which runs from -2 to -1 between
                    % minAlpha and alpha1
                    indB(k) = -1 + (alphaB(k)-alpha1)/(alpha1-minAlpha);
                    if alphaB(k) <= minAlpha
                        kB1a = k;
                    end
                    kB1b = k;
                elseif alphaB(k) < alpha2
                    % create index which runs from -1 to 0 between
                    % alpha1 and alpha2
                    indB(k) = (alphaB(k) - alpha2)/(alpha2-alpha1);
                    kB2 = k;
                elseif alphaB(k) <= alpha3
                    % create index which runs from 0 to 1 between
                    % alpha2 and alpha3
                    indB(k) = (alphaB(k) - alpha2)/(alpha3-alpha2);
                    kB3a = k+1;
                else
                    % create index which runs from 1 to 2 between
                    % alpha3 and maxAlpha
                    indB(k) = 1 + (alphaB(k)-alpha3)/(maxAlpha-alpha3);
                    if alphaB(k) <= maxAlpha
                        kB3b = k;
                    end
                end
            end
            if kB3a >= NA, kB3a=NA; kB3b=NA-1; end;
            
            ind1 = linspace(-1,0,max((kA2-kA1b),(kB2-kB1b)));
            ind2 = linspace( 0,1,max((kA3a-kA2),(kB3a-kB2)));
            ind = [indA(kA1a:kA1b);ind1';ind2(2:end)';indA(kA3a:kA3b)];
            
            alphaC = interp1(indA,alphaA,ind)*(1-w) + interp1(indB,alphaB,ind)*w;
            CLC = interp1(indA,CLA,ind)*(1-w) + interp1(indB,CLB,ind)*w;
            CDC = interp1(indA,CDA,ind)*(1-w) + interp1(indB,CDB,ind)*w;
%             LoDC = interp1(indA,CLA./CDA,ind)*(1-w) + interp1(indB,CLB./CDB,ind)*w;
%             CDC = CLC./LoDC;
            CMC = interp1(indA,CMA,ind)*(1-w) + interp1(indB,CMB,ind)*w;

            af = PolarDef();
            af.rawlist = {'alpha','CL','CD','CM'};
            af.rawdata = [alphaC,CLC,CDC,CMC];
            
            % Initialize modlist & moddata
            af.modlist = {'alpha','CL','CD','CM'};
            af.resetModData;
            af.updateIntrinsic;
        end
        
        function plotInterp(obj,varlist,X,Xq,delay,alpha1,alpha3)

            if nargin<=5
                alpha1=[];
            end
            if nargin<=6
                alpha3=[];
            end
            
            clf;
%             subplot(1,2,1);
            obj.plotMod(varlist{1},varlist{2},'k-');
            t(1) = title('');
            xlabel(varlist{1});
            ylabel(varlist{2});
            set(gca,'XLimMode','manual');
            hold on;
            hl(1) = plot(0,0,'r:');
            hp(1) = plot(0,0,'r.');
            hold off;
            
%             subplot(1,2,2);
%             obj.plotMod('alpha',varlist{2},'k-'); hold on;
%             hl(2) = plot(0,0,'r:');
%             hp(2) = plot(0,0,'r.'); hold off;
%             t(2) = title('');
%             xlabel('alpha');
%             ylabel(varlist{2});
%             set(gca,'XLimMode','manual');
            
            for Xquery=Xq
                afC = obj.interp1(X,Xquery,alpha1,alpha3);
                x = afC.getModData(varlist{1});
                y1 = afC.getModData(varlist{2});
%                 y2 = afC.getModData(varlist{2});
                dp = round(length(x)/40);
                titleStr = sprintf('Xq = %g',Xquery);
                if ishandle(hl(1))
                    set(hl(1),'XData',x,'YData',y1);
                    set(hp(1),'XData',x(1:dp:end),'YData',y1(1:dp:end));
                    set(t(1),'String',titleStr);
                end
%                 if ishandle(hl(2))
%                     set(hl(2),'XData',x,'YData',y2);
%                     set(hp(2),'XData',x(1:dp:end),'YData',y2(1:dp:end));
%                     set(t(2),'String',titleStr);
%                 end
                if (~isempty(delay) || delay>0) && all(ishandle(hl))
                    pause(delay);
                end
            end
        end
        
        function plotIntrinsic(obj)
            N = numel(obj);
            for k = 1:N
                x = obj(k).getModData('alpha');
                y = obj(k).getModData('CL');
                p = obj(k).param.intrinsic;
                if k>1, hold on; end  % "add" subsequent plots
                if isequal(k,1),holdState=ishold; end; % remember original hold state
                a = find(x>=p.alphaTrend(1),1,'first');
                b = find(x<=p.alphaTrend(2),1,'last');
                slopeInt = p.SlopeIntercept;
                h=plot(x,y,'k',...
                       x(a:b),polyval(slopeInt,x(a:b)),'m--',...
                       p.maxCLAlpha,p.maxCL,'k+',...
                       p.zeroCLAlpha,0,'k+',...
                       p.maxLoDAlpha,p.maxLoDLift,'k*');
                set(h(2),'LineWidth',2);
            end
            if ~holdState
                hold off; % reset original hold state if "off"
            end
            xlabel('Alpha');
            ylabel('CL');
        end
        
        function af_table = writePolar(obj,opts)
            af_table.title = opts.title;
            N = numel(obj);
            af_table.NumTabs = N;
            
            for k=1:N
                if isfield(obj(k).param,'source') && isfield(obj(k).param.source,'Re')
                    af_table.Tab(k).Re = obj(k).param.source.Re;
                else
                    af_table.Tab(k).Re = 9.99e9;
                end
                if isfield(obj(k).param,'source') && isfield(obj(k).param.source,'Ctrl')
                    af_table.Tab(k).Ctrl = obj(k).param.source.Ctrl;
                else
                    af_table.Tab(k).Ctrl = 0;
                end
                if isfield(obj(k).param,'AeroDyn')
                    af_table.Tab(k).AlfaStal = obj(k).param.AeroDyn.AlfaStal;
                    af_table.Tab(k).AOL      = obj(k).param.AeroDyn.AOL;
                    af_table.Tab(k).CnA      = obj(k).param.AeroDyn.CnA;
                    af_table.Tab(k).CnS      = obj(k).param.AeroDyn.CnS;
                    af_table.Tab(k).CnSL     = obj(k).param.AeroDyn.CnSL;
                    af_table.Tab(k).AOD      = obj(k).param.AeroDyn.AOD;
                    af_table.Tab(k).Cd0      = obj(k).param.AeroDyn.Cd0;
                else
                    error('Need DynStall parameters to write AeroDyn polar.');
                end
                af_table.Tab(k).NumAlf = size(obj(k).moddata,1);
                af_table.Tab(k).Alpha  = obj(k).getModData('alpha');
                af_table.Tab(k).Cl     = obj(k).getModData('CL');
                af_table.Tab(k).Cd     = obj(k).getModData('CD');
                if isfield(opts,'UseCM') && opts.UseCM
                    af_table.Tab(k).Cm     = obj(k).getModData('CM');
                end
            end
        end
        
    end
    
    methods (Static)
        function obj = cylinderPolar(varargin)
            obj = PolarDef();
            obj.rawlist = {'alpha','CL','CD','CM'};
            alpha = transpose(-180:10:180);
            CL =  0.0001*ones(size(alpha));
            CD =  0.3500*ones(size(alpha));
            CM = -0.0001*ones(size(alpha));
            obj.rawdata = [alpha,CL,CD,CM];
            obj.modlist = {'alpha','CL','CD','CM'};
            obj.resetModData;
            obj.param.intrinsic.manual = true;
            for k=2:2:nargin
                label = varargin{k-1};
                value = varargin{k};
                if ischar(label) && isfield(obj.param.intrinsic,label)
                    if isscalar(value) && isfloat(value)
                        obj.param.intrinsic.(label) = value;
                    end
                end
            end
            
        end
        
        function mopt = createModOpts(varargin)
            pairedArgs = isequal(0,rem(nargin,2));
            assert(pairedArgs,'PolarDef:nargchk','createModOpt requires PAIRED key,value arguments');
            
            charKeys = cellfun(@(x) ischar(x), varargin(1:2:end));
            assert(all(charKeys),'PolarDef:argtype','createModOpt requires the key in a key,value pair to be a STRING');
            
            mopt = struct(varargin{:});
        end
        
        function zeroCLAlpha = findZeroCL(alpha,CL)
            N = length(CL);
            k = 1;
            while(k<N && alpha(k)<0)
                k=k+1; % increase to where alpha==0
            end
            while (k>0 && CL(k)>0)
                k=k-1; % decrease to where CL==0
            end
            k=k+1; % want k where CL(k-1) < 0 < CL(k)
            if k==1
                if N>=2
                    % approximate with slope from last two points
                    slope = (CL(2)-CL(1))/(alpha(2)-alpha(1));
                    zeroCLAlpha = alpha(1)-CL(1)/slope;
                else
                    zeroCLAlpha = alpha(1);
                end
            else
                w = CL(k-1)/(CL(k-1)-CL(k));
                zeroCLAlpha = alpha(k-1)*(1-w) + alpha(k)*w;
            end
        end
        
        function [maxCL,maxCLAlpha] = findMaxCL(alpha,CL)
            N = length(alpha);
            k = 1;
            while(k<N && alpha(k)<0)
                k=k+1;
            end
            maxCL = CL(k);
            maxCLAlpha = alpha(k);
            k=k+1;
            while(k<=N)
                if CL(k)>=maxCL
                    maxCL = CL(k);
                    maxCLAlpha = alpha(k);
                else
                    % CL has decreased
                    break; % use only the first local max found
                end
                k=k+1;
            end
        end
        
        function [maxLoD,maxLoDAlpha] = findMaxLoD(alpha,LOD)
            [maxLoD,k] = max(LOD);
            maxLoDAlpha = alpha(k);
        end
        
        function [alphaTrend,SlopeIntercept] = findAlphaTrend(alpha,CL)
            alpha = alpha(:);
            CL = CL(:);
            N = length(alpha);
            k = 1;
            while(k<N && alpha(k)<0)
                k=k+1;
            end
            diffSlope = diff(CL) ./ diff(alpha);
            tol = 0.05;
            slopeLB = diffSlope(k)*(1-tol);
            slopeUB = diffSlope(k)*(1+tol);
            if k>1
                a = k-1;
                b = k;
            else
                a = k;
                b = k;
            end
            while (a>1 && diffSlope(a)>slopeLB && diffSlope(a)<slopeUB)
                a=a-1;
            end
            while (b<N-1 && diffSlope(b)>slopeLB && diffSlope(b)<slopeUB)
                b=b+1;
            end
            alphaTrend = [alpha(a) alpha(b)];
%             H = [alpha(a:b), ones(b-a+1,1)];
%             P = (H'*H)\H'*CL(a:b);
            P = polyfit(alpha(a:b),CL(a:b),1);
            SlopeIntercept = P;
        end
        
        function recolorplot(StackIndex,ColorOrder)
            % ``recolorplot(StackIndex,ColorOrder)``
            % where ``StackIndex`` specifies which axes children to change 
            % (Note that the last child to be added shows up first on the "stack")
            % ``ColorOrder`` is the list of colors to use
            % Modify last six axes children (top six on stack):
            %
            % Example:
            %
            %	``recolorplot(1:6,jet(6));``
            %
            h = get(gca,'Children');
            assert(all(StackIndex>0),'all N must be positive');
            assert(max(StackIndex)<=numel(h),'max(N) must be <= number of axes children');
            N = size(ColorOrder,1);
            for k=numel(StackIndex):-1:1
                kmod = mod(StackIndex(k)-1,N)+1;
                set(h(k),'Color',ColorOrder(kmod,:));
            end
        end
    end
    
    methods (Hidden)
        function readXFoil(obj,filelines)
            % read Xfoil polar (works with v6.96)
            m = regexp(filelines{2},'\S+','match');
            obj.source{1} = m{1};  % XFOIL
            obj.source{2} = m{3};  % version number
            obj.titleLine = strtrim(filelines{4});  % "Calculated for..."
            obj.notes{1}  = strtrim(filelines{6});
            m = regexp(filelines{8},'[\d.-+]+','match');
            obj.param.source.xtrf_top    = str2double(m{1});
            obj.param.source.xtrf_bottom = str2double(m{2});
            m = regexp(filelines{9},'[\d.-+]+','match');
            obj.param.source.Mach  = str2double(m{1});
            obj.param.source.Re    = str2double(m{2}) * 10^str2double(m{3});
            obj.param.source.Ncrit = str2double(m{4});
            obj.rawlist = regexp(filelines{11},'\S+','match');
            nheaderlines = 12;
            nrows = numel(filelines) - nheaderlines;
            ncols = numel(obj.rawlist);
            obj.rawdata = zeros(nrows,ncols);
            for k = 1:nrows
                d = sscanf(filelines{nheaderlines+k},'%f');
                if isequal(ncols,numel(d));
                    obj.rawdata(k,:) = transpose(d);
                else
                    obj.rawdata(k:end,:) = [];
                    break;
                end
            end
        end
        
        function readCactus(obj,filelines)
            % read Cactus polar data (first table only)
            obj.titleLine = strtrim(filelines{1});
            obj.source{1} = 'CACTUS';
            obj.notes{1}  = strtrim(filelines{2});
            obj.notes{2}  = strtrim(filelines{3});
            obj.notes{3}  = strtrim(filelines{4});
            
            Nlines = numel(filelines);
            maxTables = 1;
            tableRows = zeros(maxTables,2);
            inTable = false;
            j = 1; 
            for k=1:Nlines
                if ~inTable
                    if strfind(filelines{k},'Reynolds')
                        tableRows(j,1)=k;  % first line (Reynolds Number)
                        inTable = true;
                    end
                else
                    if length(filelines{k})==1
                        tableRows(j,2)=k-1; % last line
                        inTable = false;
                        j = j + 1;
                    end
                end
                if j>maxTables
                    break;
                end
            end
            
            tableFirstLine = tableRows(1,1);
            n = regexp(filelines{tableFirstLine},'.*:(?<Re>.*)','names');
            obj.param.source.Re = str2double(strtrim(n(1).Re));
            obj.rawlist = {'alpha','CL','CD','CM'};
            
            nheaderlines = 7;
            nrows = tableRows(1,2) - tableRows(1,1) - nheaderlines + 1;
            ncols = numel(obj.rawlist);
            obj.rawdata = zeros(nrows,ncols);
            offset = tableRows(1,1)+nheaderlines-1;
            for k = 1:nrows
                d = sscanf(filelines{offset+k},'%f');
                if isequal(ncols,numel(d));
                    obj.rawdata(k,:) = transpose(d);
                else
                    obj.rawdata(k:end,:) = [];
                    break;
                end
            end
        end
        
        function readAerodyn(obj,filelines)
            % read Aerodyn polar data (first table only)
            aerodynVersion = regexp(filelines{1},'v\d{2}\.*\d*','match','once');
            switch aerodynVersion
                case {'v13','v13.0','13.0'}
                    ADver = 'v13';
                case {'v12','v12.0','12.0'}
                    ADver = 'v12';
                otherwise
                    ADver = 'v12';
            end
            obj.titleLine = '';
            obj.source{1} = 'AERODYN';
            
            switch ADver
                case 'v13'
                    obj.notes{1}  = strtrim(filelines{2});
                    obj.notes{2}  = strtrim(filelines{3});
%                     nTables = sscanf(filelines{4},'%g',[1 1]);
                    obj.param.source.Re = 1e6 * sscanf(filelines{5},'%g',[1 1]);  % Reynolds number in millions
                    obj.param.source.Ctrl = sscanf(filelines{6},'%g',[1 1]);
                    obj.param.source.AlfaStal = sscanf(filelines{7},'%g',[1 1]);
                    obj.param.source.AOL  = sscanf(filelines{8},'%g',[1 1]);
                    obj.param.source.CnA  = sscanf(filelines{9},'%g',[1 1]);
                    obj.param.source.CnS  = sscanf(filelines{10},'%g',[1 1]);
                    obj.param.source.CnSL = sscanf(filelines{11},'%g',[1 1]);
                    obj.param.source.AOD  = sscanf(filelines{12},'%g',[1 1]);
                    obj.param.source.Cd0  = sscanf(filelines{13},'%g',[1 1]);
                    offset = 13;
                case 'v12'
                    obj.notes{1}  = strtrim(filelines{1});
                    obj.notes{2}  = strtrim(filelines{2});
                    nTables = sscanf(filelines{3},'%g',[1 1]);
                    obj.param.source.TableID = sscanf(filelines{4},'%g',[1 1]);
                    obj.param.source.AlfaStal = sscanf(filelines{5},'%g',[1 1]);
                    obj.param.source.AOL  = sscanf(filelines{9},'%g',[1 1]);
                    obj.param.source.CnA  = sscanf(filelines{10},'%g',[1 1]);
                    obj.param.source.CnS  = sscanf(filelines{11},'%g',[1 1]);
                    obj.param.source.CnSL = sscanf(filelines{12},'%g',[1 1]);
                    obj.param.source.AOD  = sscanf(filelines{13},'%g',[1 1]);
                    obj.param.source.Cd0  = sscanf(filelines{14},'%g',[1 1]);
                    offset = 14; 
            end
            
            k=1;
            while true
                d = sscanf(filelines{offset+k},'%f');
                if isempty(d)
                    break;
                end
                if k==1
                    dataArray = transpose(d);
                else
                    dataArray(k,:) = transpose(d);
                end
                k=k+1;
            end

            nCols = size(dataArray,2);
            switch ADver
                case 'v13'
                    obj.rawdata = dataArray;
                    if nCols==4
                        obj.rawlist = {'alpha','CL','CD','CM'};
                    elseif nCols==3
                        obj.rawlist = {'alpha','CL','CD'};
                    else
                        error('Expecting 3 or 4 column table');
                    end 
                case 'v12'
                    if (nCols-1)/nTables==3
                        obj.rawlist = {'alpha','CL','CD','CM'};
                        obj.rawdata = dataArray(:,1:4);
                    elseif (nCols-1)/nTables==2
                        obj.rawlist = {'alpha','CL','CD'};
                        obj.rawdata = dataArray(:,1:3);
                    else
                        error('Unexpected number of columns');
                    end
            end
            
        end
        
    end
    
    % END OF CLASSDEF
end

function fcontents = freadcontents(filename)
    fid = fopen(filename);
    assert(fid ~= -1,'FileIO:FileOpenError',...
        'Could not open file "%s"',filename);
    fcontents = fread(fid,inf,'*char')';
    fclose(fid);
end

%==========================================================================

function [Table3D,Trend] = Selig_3DStall(Table2D,RPM,R,V,rOverR,Chord,AlphaMinTrend,AlphaMaxTrend,AlphaEnd,abd)
% Program to apply Du-Selig and Eggars 3-D stall corrections to 2-D airfoil data
%   C Hansen, Windward Engineering, Dec 2003
%   Reproduced from AirfoilPrep_v2p2.xls available at wind.nrel.gov
%
% Table3D = stall3d(Table2D,RPM,R,V,rOverR,Chord,AlphaMinTrend,AlphaMaxTrend,AlphaEnd)

if ~exist('abd','var') || isempty(abd)
    a=1; b=1; d=1;  % Basic constants of Selig method.  Do not change except for research
else
    a=abd(1);
    b=abd(2);
    d=abd(3);
end

DToR = pi/180;  % DToR still needed for FL calculation

% "Call Update2DSlope", or rather perform equivalent steps
nTable2D = size(Table2D,1);
Alpha2D = Table2D(:,1);  % first column Alpha
CL2D = Table2D(:,2);     % second column CL
CD2D = Table2D(:,3);     % third column CD
rows = ((Alpha2D > AlphaMinTrend) & (Alpha2D < AlphaMaxTrend));
TrendInput = Table2D(rows,1:2); % select rows between AlphaMin and AlphaMax
% perform linear fit to calculate slope and intercept
P = polyfit(TrendInput(:,1),TrendInput(:,2),1); % linear fit
% calculate some constants, uses original notation for Selig variables
CLSlope = P(1);     % get slope
CLIntercept = P(2); % get intercept
AlphaZero = -CLIntercept/CLSlope;
Trend = [TrendInput(:,1),polyval(P,TrendInput(:,1))];

% calculate some constants
cOverr = Chord / R / rOverR;  % c/r
Lambda = RPM * pi/30 * R / sqrt(V^2 + (RPM * pi/30 * R)^2);
Expon = d / Lambda / rOverR;
FL = 1 / (CLSlope / DToR) * ((1.6 * cOverr / 0.1267)*(a - cOverr ^ Expon) / (b + cOverr ^ Expon) - 1);

% calculate 3D values
AlphaRad = DToR * Alpha2D;
%CN2D = CL2D .* cos(AlphaRad) + CD2D .* sin(AlphaRad);
%CT2D = CL2D .* sin(AlphaRad) - CD2D .* cos(AlphaRad);

CLP = CLSlope * (Alpha2D - AlphaZero);
DelCL1 = FL * (CLP - CL2D);

adj = ones(nTable2D,1);
rows = (Alpha2D > AlphaEnd);
adj(rows) = ((90 - Alpha2D(rows)) ./ (90 - AlphaEnd)).^2;

CL3D = CL2D + DelCL1 .* adj;
DelCL2 = CL3D - CL2D;

DelCD = DelCL2 .* (sin(AlphaRad) - 0.12 * cos(AlphaRad)) ./ (cos(AlphaRad) + 0.12 * sin(AlphaRad));
CD3D = CD2D + DelCD;

% Return the results if in the appropriate range
%  (note this method does not work for very large or small angles)
Table3D = [Alpha2D, CL3D, CD3D];
rows = ((Alpha2D >= -10) & (Alpha2D <= 90));
Table3D = Table3D(rows,:);

end


%==========================================================================
function OutputTable = TableExtrap(InputTable,CDMax,UseCM,AlphaFP)
%Foilcheck program converted to Visual Basic routine
% CH Windward Engineering, Dec 2003
%
% OutputTable = TableExtrap(InputTable,CDMax,UseCM,[AlphaFP])

if ~exist('AlphaFP','var') || isempty(AlphaFP)
    UseFlatPlate = false;
else
    % get CD from flat plate theory for abs(angles) > AlphaFP
    UseFlatPlate = true;
end

DToR = pi/180;

% read in original table values
nTable1 = size(InputTable,1);
Alpha1 = InputTable(:,1);
CL1 = InputTable(:,2);
CD1 = InputTable(:,3);
if (UseCM)
    CM1 = InputTable(:,4);
else
    CM1 = zeros(size(Alpha1));
end

% apply viterna extrapolation
% ==============================
%   CLRef        % CL at upper matching point
%   CMCoef       % Coefficient used in CM calcs
%   VAlphaLo     % lower matching point
%   VAlphaHi     % upper matching point
VAlphaLo = Alpha1(1);     % use first and last points for matching!
VAlphaHi = Alpha1(end);

% Find CL and CD at VAlphaHi, the upper matching point
for i = 1:nTable1
    if (abs(Alpha1(i) - VAlphaHi) < 0.01)
        CLRef = CL1(i);
        CDRef = CD1(i);
        %fprintf(' CLRef = %g\n',CLRef);
        break
    end
    if (i == nTable1)  % we shouldn't get here
        error('Error: Your upper matching point for the Viterna calculation was not found in the original airfoil table');
    end
end

% get Viterna coefficients
SAlpha = sin(VAlphaHi * DToR);
CAlpha = cos(VAlphaHi * DToR);

A2 = (CLRef - CDMax / 2 * sin(2 * VAlphaHi * DToR)) * SAlpha / CAlpha ^ 2;
B2 = (CDRef - CDMax * SAlpha ^ 2) / CAlpha;

% get pitching moment coefficients (if needed)
if (UseCM)
    % call CMCoeff
    hi = find(Alpha1 >= VAlphaHi,1); % jcb: this should be same as numel(Alpha1)
    CLHi = CL1(hi);  % jcb: or CL1(end)
    CDHi = CD1(hi);  % jcb: or CD1(end)
    CMHi = CM1(hi);  % jcb: or CM1(end)
    
    FoundZeroLift = false;
    % get CM at angle of zero lift (CM0)
    for i = 1:nTable1-1
        if (abs(Alpha1(i)) < 20 && CL1(i) <= 0 && CL1(i+1) >= 0)
            p = -CL1(i) / (CL1(i+1) - CL1(i));
            CM0 = CM1(i) + p*(CM1(i+1) - CM1(i));
            FoundZeroLift = true;
            break
        end
    end
    if (~FoundZeroLift)  % zero lift not in range of orig table, use first two points
        p = -CL1(1) / (CL1(2) - CL1(1));  % interpolation parameter
        CM0 = CM1(1) + p*(CM1(2) - CM1(1));
    end
    
    XM = (-CMHi + CM0) / (CLHi * cos(VAlphaHi * DToR) + CDHi * sin(VAlphaHi * DToR));
    CMCoef = (XM - 0.25) / tan((VAlphaHi - 90) * DToR);
    
end

% Create new table that cover entire AoA range. Fill in CL, CD, CM in original range
% call NewTable
Step = 10;   % angle of attack step in deg
SizeMax = 500;  % dimension of Alpha arrays
HiEnd = 0;      % flag says we've gone past VAlphaHi
remainder = VAlphaHi - floor(VAlphaHi / Step) * Step;

Alpha2 = zeros(SizeMax,1);
CL2 = zeros(SizeMax,1);
CD2 = zeros(SizeMax,1);
CM2 = zeros(SizeMax,1);

Alpha2(1) = -180;
count = 1;
nTable2 = SizeMax;
for i = 2:SizeMax
    Alpha2(i) = Alpha2(i-1) + Step;
    if (Alpha2(i) > 180 - Step)         % stop when we reach 180 deg
        Alpha2(i) = 180;
        nTable2 = i;
        break
    end
    if (Alpha2(i) >= VAlphaLo)
        if (Alpha2(i) - Step < VAlphaHi) %in the range of the original table, use those alphas
            Alpha2(i) = Alpha1(count);    %save these values and write to worksheet
            CL2(i) = CL1(count);
            CD2(i) = CD1(count);
            CM2(i) = CM1(count);
            count = count + 1;    % number of original alpha values we've used
        elseif (HiEnd == 0)   % we just left the upper end of the original alpha table
            HiEnd = 1;
            Alpha2(i) = VAlphaHi - remainder + Step; %find first multiple of 10 after ValphaHi
        end
    end
end

% Now fill new table with calculated values that are outside original table range
% call ViternaFill
CLAdj = 0.7;  % adjustment factor for CL when abs(Alpha)>90

lo = find(Alpha1 >= VAlphaLo,1);
hi = find(Alpha1 >= VAlphaHi,1);
CLLo = InputTable(lo,2);  % values at each end of original table
CDLo = InputTable(lo,3);
CLHi = InputTable(hi,2);
CDHi = InputTable(hi,3);

%Get angle when flap plate theory will be applied (if that option selected)
%%AlphaFP = Cells(11, 2).Value

for i = 1:nTable2    %step through all alphas.  Notice that if angle is
    %within original table range the value is not changed
    Alfa = Alpha2(i);
    SAlfa = sin(Alfa * DToR);
    CAlfa = cos(Alfa * DToR);
    if (Alfa > 180 || Alfa < -180)
        error(['Angle of attack = ',Alfa,' outside range + to -180 deg in Viterna calculation'])
    elseif (Alfa >= VAlphaHi && Alfa <= 90)
        CL2(i) = CDMax / 2 * sin(2 * Alfa * DToR) + A2 * CAlfa ^ 2 / SAlfa;
        CD2(i) = CDMax * SAlfa ^ 2 + B2 * CAlfa;
    elseif (Alfa > 90 && Alfa <= 180 - VAlphaHi)
        Ang = 180 - Alfa;
        Sang = sin(Ang * DToR);
        Cang = cos(Ang * DToR);
        CL2(i) = CLAdj * (-CDMax / 2 * sin(2 * Ang * DToR) - A2 * Cang ^ 2 / Sang);
        CD2(i) = CDMax * Sang ^ 2 + B2 * Cang;
    elseif (Alfa > 180 - VAlphaHi && Alfa <= 180)
        Ang = Alfa - 180;
        Sang = sin(Ang * DToR);
        Cang = cos(Ang * DToR);
        CL2(i) = Ang / VAlphaHi * CLHi * CLAdj;
        CD2(i) = CDMax * Sang ^ 2 + B2 * Cang;
    elseif (Alfa > -VAlphaHi && Alfa < VAlphaLo)
%         CL2(i) = CLAdj * (-CLHi + (Alfa + VAlphaHi) / (VAlphaHi + VAlphaLo) * (CLHi + CLLo));
%jcb 2014-03-12: the original equation (above) results in 
%                CL2(VAlphaLo) = CLAdj*CLLo   rather than CLLo 
        CL2(i) = -CLHi*CLAdj + (Alfa + VAlphaHi) / (VAlphaHi + VAlphaLo) * (CLHi*CLAdj + CLLo);
        CD2(i) = CDLo + (-Alfa + VAlphaLo) * (CDHi - CDLo) / (VAlphaHi + VAlphaLo);
    elseif (Alfa <= -VAlphaHi && Alfa >= -90)
        Ang = -Alfa;
        Sang = sin(Ang * DToR);
        Cang = cos(Ang * DToR);
        CL2(i) = CLAdj * (-CDMax / 2 * sin(2 * Ang * DToR) - A2 * Cang ^ 2 / Sang);
        CD2(i) = CDMax * Sang ^ 2 + B2 * Cang;
    elseif (Alfa < -90 && Alfa >= -180 + VAlphaHi)
        Ang = 180 + Alfa;
        Sang = sin(Ang * DToR);
        Cang = cos(Ang * DToR);
        CL2(i) = CLAdj * (CDMax / 2 * sin(2 * Ang * DToR) + A2 * Cang ^ 2 / Sang);
        CD2(i) = CDMax * Sang ^ 2 + B2 * Cang;
    elseif (Alfa < -180 + VAlphaHi && Alfa >= -180)
        Ang = 180 + Alfa;
        Sang = sin(Ang * DToR);
        Cang = cos(Ang * DToR);
        CL2(i) = Ang / VAlphaHi * CLHi * CLAdj;
        CD2(i) = CDMax * Sang ^ 2 + B2 * Cang;
    end % if
    
    if (CD2(i) < 0)   %watch out for negative CD's
        CD2(i) = 0.01;
    end
end

%Change to using CD from flat plate theory if desired
if (UseFlatPlate)
    for i = 1:nTable2
        Alfa = Alpha2(i);
        if (abs(Alfa) >= AlphaFP)
            CD2(i) = CL2(i) * tan(Alfa * DToR);
        end
    end
end

%Do CM calculations and write value to output table
if (UseCM)
    CMLo = CM1(1);
    for i = 1:nTable2
        Alfa = Alpha2(i);
        if (Alfa >= VAlphaLo && Alfa <= VAlphaHi)
            continue %no action needed
        end
        if (Alfa > -165 && Alfa < 165)
            if (abs(Alfa) < 0.01)
                CM2(i) = CM0;
            else
%                 if (Alfa > 0)
%                     x = CMCoef * tan((Alfa - 90) * DToR) + 0.25; % coefficient used in CM calc
%                     CM2(i) =   CM0 - x * (CL2(i) * cos(Alfa * DToR) + CD2(i) * sin(Alfa * DToR));
%                 else
%                     % correction by Scott Larwood for negative Alfa (AirfoilPrep_v2.2.1)
%                     x = CMCoef * tan((-Alfa - 90) * DToR) + 0.25; % coefficient used in CM calc
%                     CM2(i) = -(CM0 - x * (-CL2(i) * cos(-Alfa * DToR) + CD2(i) * sin(-Alfa * DToR)));
%                 end

                %jcb: this creates a continuous curve at VAlphaLo, but it
                % is just a shot in the dark. I have no idea whether it is
                % correct.
                if (Alfa > 0)
                    XM = (-CMHi + CM0) / (CLHi * cos(VAlphaHi * DToR) + CDHi * sin(VAlphaHi * DToR));
                    CMCoef = (XM - 0.25) / tan((VAlphaHi - 90) * DToR);
                else
                    XM = (-CMLo + CM0) / (CLLo * cos(VAlphaLo * DToR) + CDLo * sin(VAlphaLo * DToR));
                    CMCoef = (XM - 0.25) / tan((VAlphaLo - 90) * DToR);
                end
                x = CMCoef * tan((Alfa - 90) * DToR) + 0.25; % coefficient used in CM calc
                CM2(i) =   CM0 - x * (CL2(i) * cos(Alfa * DToR) + CD2(i) * sin(Alfa * DToR));
            end
        else
            switch (Alfa)
                case 165
                    CM2(i) = -0.4;
                case 170
                    CM2(i) = -0.5;
                case 175
                    CM2(i) = -0.25;
                case 180
                    CM2(i) = 0;
                case -165
                    CM2(i) = 0.35;
                case -170
                    CM2(i) = 0.4;
                case -175
                    CM2(i) = 0.2;
                case -180
                    CM2(i) = 0;
                otherwise
                    error('Angle encountered for which there is no CM table value (near +/-180°).  Program will stop');
            end
        end
    end %for i
    
    OutputTable = [Alpha2(1:nTable2), CL2(1:nTable2), CD2(1:nTable2), CM2(1:nTable2)];
    
else
    
    OutputTable = [Alpha2(1:nTable2), CL2(1:nTable2), CD2(1:nTable2)];
    
end %if (UseCM)

end % function


%==========================================================================
function params = DynStall(InputTable,AlphaMinTrend,AlphaMaxTrend,StallAngle,NegStallCn)
%Foilcheck program converted to Visual Basic routine
% CH Windward Engineering, Dec 2003
%
% StallHeader = DynStall(InputTable,AlphaMinTrend,AlphaMaxTrend,StallAngle,NegStallCn)

DToR = pi/180;

Alpha1 = InputTable(:,1);
CL1 = InputTable(:,2);
CD1 = InputTable(:,3);
% CM1 = InputTable(:,4);
CN = CL1 .* cos(DToR * Alpha1) + CD1 .* sin(DToR * Alpha1);
% CT = CL1 .* sin(DToR * Alpha1) - CD1 .* cos(DToR * Alpha1);

% get CN slope
rows = ((Alpha1 > AlphaMinTrend) & (Alpha1 < AlphaMaxTrend));
CNSlopeTable = [Alpha1(rows),CN(rows)];
P = polyfit(CNSlopeTable(:,1),CNSlopeTable(:,2),1);
CNSlope = P(1)*180/pi;           % Cn slope for zero lift
ZeroCN = -P(2)/CNSlope*180/pi;   % Zero Cn angle of attack (deg)
CNExtrap = CNSlope*pi/180*(StallAngle-ZeroCN);

% get angle of attack for CDmin (search only in range -20 to +20 deg)
rows = find(((Alpha1 > -20) & (Alpha1 < 20)));
[~,i] = min(CD1(rows));
CDMin = CD1(rows(i));
AlphaCDMin = Alpha1(rows(i));

params.AlfaStal = StallAngle;
params.AOL  = ZeroCN;
params.CnA  = CNSlope;
params.CnS  = CNExtrap;
params.CnSL = NegStallCn;
params.AOD  = AlphaCDMin;
params.Cd0  = CDMin;

end