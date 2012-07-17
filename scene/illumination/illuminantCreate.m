function lightSource = illuminantCreate(iName,temperature,luminance,spectrum)
% Create a light source structure.  Only used for Macbeth chart creation
%
%  lightSource = illuminantCreate(iName,temperature,luminance,spectrum)
%
% Create a lightSource data structure, containing information about an
% illuminant.  The lightSource is organized approximately as a scene
% structure and can be addressed using sceneGet.  But it has a unique
% status.  I don't think it is used anywhere else, so it has not been
% thoroughly tested.
%
% A few standard light sources are supported at present.  These are d65,
% d50, tungsten, fluorescene, blackbody's (over a range of color
% temperatures), 550nm.  See the internal routine, readIllumination (below). 
%
% Examples:
%   lightSource = illuminantCreate('d65')
%   lightSource = illuminantCreate('blackbody',3500,100)
%   lightSource = illuminantCreate('blackbody',6500,100)
%
%   spectrum.wave = (380:4:1068);
%   lightSource = illuminantCreate('equalEnergy',[],100,spectrum)
%
% Copyright ImagEval Consultants, LLC, 2005.

% TODO
% This code should be updated.  It is only used in creating the Macbeth
% Chart at present.  The assumption made here is that the data on disk are
% in energy units. The data are then converted to photons.  In the future,
% the units of the file will be stored.
if ieNotDefined('iName'),       error('You must specify a light source name.'); end
if ieNotDefined('temperature'), temperature = []; end
if ieNotDefined('luminance'),   luminance = 100; end

lightSource.name = iName;
lightSource.type = 'scene';

lightSource.temperature = temperature;	        % [K]
lightSource.luminance   = luminance;		    % [cd/m^2]
if ieNotDefined('spectrum'), lightSource = initDefaultSpectrum(lightSource,'hyperspectral');
else                         lightSource.spectrum = spectrum;
end

wave = sceneGet(lightSource,'wave');
illuminantEnergy = readIllumination(lightSource);		% [W/(sr m^2 nm)]
illuminantPhotons = Energy2Quanta(wave,illuminantEnergy); % Check this step

if strcmp('blackbody',iName)
    lightSource.name = sprintf('blackbody%.0f',lightSource.temperature);
end

% Creating more of the lightSource structure
lightSource.data.photons = illuminantPhotons;% [photons/(s sr m^2 nm)]

return;
