function Out = nrlmsise2_1(lat, lon, time, alt, varargin)
% MATLAB wrapper for the NRLMSISE2.1 model and MEX function
%
% INPUTS
%   latitude(deg), longitude(deg), time(datenum or datetime), altitude(km)
%   latitude and longitude must have the same size, unless the 'map'
%   keyword if provided. 
%
%   time must be a scalar.
%
% OUTPUTS
%   structure with fields
%   altitude(km), latitude(deg), longitude(deg),
%   electron density (Ne/m^3) and F107 flux
%       TN:    Temperature at altitude (K)
%       Total: Total mass density (kg/m3)
%       N2:    N2 number density (m-3)
%       O2:    O2 number density (m-3)
%       O:     O number density (m-3)
%       He:    He number density (m-3)
%       H:     H number density (m-3)
%       Ar:    Ar number density (m-3)
%       N:     N number density (m-3)
%       Oa:    Anomalous oxygen number density (m-3)
%       NO:    NO number density (m-3)
%       TEx:   Exospheric temperature (K) (optional argument)
%
% OPTIONAL INPUTS
%
%   By default, NRLMSISE works in profile mode. Lat and Lon must be the
%   same size, and the output will calculate NRLMSISE2.1 at every altitude 
%   for each Lat/Lon pair. 
%
%   'map' mode: in this mode, all unique values of inputs are expanded in a
%   4d grid. The output Ne is of a size [n lat,n lon,n alt,n time]
%
%   NRLMSISE2.1 will read F107 flux from a file. Use the 'flux'
%   keyword to provide your own
%   e.g. nrlmsise2_1(45,-78,...,'flux',131.2)
%
%   use the 'update' keyword to update apf107 automatically
%
%   Ben Reid 2020, 2022

ErrorStr = 'MATLAB:NRLMSISE2_1';

[ModelDir, ~, ~] = fileparts(mfilename('fullpath'));
ModelDir = fullfile(ModelDir);

if isempty(ModelDir)
    error(ErrorStr,[...
        'NRLMSISE model directory not configured! Please edit',...
        mfilename('fullpath'),' to specify the appropriate path'])
elseif ~exist(ModelDir,'dir')
    error(ErrorStr,[...
        'Invalid model directory. Directory does not exist. Please edit',...
        mfilename('fullpath'),' to specify the appropriate path'])
end

try
    OldDir = pwd;
    cd(ModelDir);
    if ~exist('mex_nrlmsise2_1','file') 

        mex -compatibleArrayDims -DO3 -Dstatic -DfPIC -Dffast-math ...
            -Dmarch=native ...
            -output mex_nrlmsise2_1 ...
            msis_constants.F90 msis_utils.F90 msis_init.F90 msis_gfn.F90 ...
            msis_tfn.F90 msis_dfn.F90 msis_calc.F90 mex_nrlmsise2_1.F
    end

    if numel(time) >1
        error(ErrorStr,[...
            'time must be a scalar']);
    end

    mode='default';
    flux = 'file';
    Update = false;
    i=1;
    while i <= numel(varargin)
        if ischar(varargin{i}) || isstring(varargin{i})
            if strcmpi(varargin{i},'default')
                mode = default;
            elseif strcmpi(varargin{i},'map')
                mode='map';
            elseif strcmpi(varargin{i},'flux')
                flux = 'input';
                i=i+1;
                F107 = varargin{i};
            elseif strcmpi(varargin{i},'update')
                Update=true;
            else
                error(ErrorStr,[...
                    'Invalid option: ',varargin{i}]);
            end
        else
            error(ErrorStr,[...
                'Invalid option.'])
        end
        i=i+1;
    end

    if Update % update apf107.dat file in iri directory
        read_apf(fullfile(ModelDir,'apf107.dat'),true);
    end

    if any((lat(:)<-90) | (lat(:)>90))
        error(ErrorStr,[...
            'Input latitude out of range [-90,90]'])
    end

    if any((lon(:)<-180) | (lon(:)>360))
        error(ErrorStr,[...
            'Input longitude out of range [-180,360]'])
    end

    APF = read_apf(fullfile(ModelDir,'apf107.dat'));

    idate = sum(datenum(APF.date(:)')<datenum(time(:))',2);

    if strncmpi(mode,'default',3)

        tvec=[numel(lat),numel(lon)];
        tvec=tvec==tvec';
        if ~all(tvec(:))
            error('In "default" mode and lat, lon  must all have the same size')
        end

        Out.time = time(:);
        Out.lon = lon;
        Out.lat = lat;
        Out.alt = alt(:);

        OutFmt = @(X) reshape(X,[numel(alt),size(lon)]);
        [~,lat] = ndgrid(alt,lat);
        [alt, lon] = ndgrid(alt,lon);

    elseif strncmpi(mode,'map',3)

        lat = unique(lat);
        lon = unique(lon);
        alt = unique(alt);

        [Out.lat,Out.lon] = ndgrid(lat,lon);
        Out.alt = alt(:);
        Out.time = time;

        OutFmt = @(X) reshape(X,numel(alt), numel(lat), numel(lon));
        [alt,lat,lon] = ndgrid(alt,lat(:),lon(:));

    end
    [Year,~,~,~,~,~] = datevec(time(:));

    IYD = floor(datenum(time(:))-datenum(Year,0,0));
    SEC = mod(datenum(time(:)),1)*24*60*60;
    F107A = APF.F107_81(idate(:));
    
    % overwrite f107from file
    if isnumeric(flux)
   
    else
        F107 = APF.F107(idate(:));
    end

    Ap = APF.Ap(idate(:));
    SW = ones(25,1);


    [D,T] = mex_nrlmsise2_1( ...
        IYD,SEC,F107A,F107,Ap,SW, ...
        [reshape(alt,1,[]); ...
        reshape(lat,1,[]); ...
        reshape(lon,1,[])]);

    Out.Total = OutFmt(D(1,:));
    Out.N2 = OutFmt(D(2,:));
    Out.O2 = OutFmt(D(3,:));
    Out.O = OutFmt(D(4,:));
    Out.He = OutFmt(D(5,:));
    Out.H = OutFmt(D(6,:));
    Out.Ar = OutFmt(D(7,:));
    Out.N = OutFmt(D(8,:));
    Out.Oa = OutFmt(D(9,:));
    Out.NO = OutFmt(D(10,:));

    Out.TEx = OutFmt(T(2,:));
    Out.TN = OutFmt(T(1,:));
    Out.F107 = F107;
    Out.F107A = F107A;
    Out.Ap = Ap;

    cd(OldDir)

catch err
    cd(OldDir)
    rethrow(err)

end




