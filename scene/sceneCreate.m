function [scene,parms] = sceneCreate(sceneName,varargin)
% Create a scene structure.
%
%  [scene,parms] = sceneCreate(sceneName,varargin)
%
% A scene describes the light emitted from a planar object, say a flat
% screen display.  The scene is located at some distance from the center of
% the optics, has a field of view, and a spectral radiance distribution.
% All of these properties, and several others, can be set. A variety of
% images can be created automatically.  The routines that create these
% scenes, including this one, serve as a template for creating others.
%
% MACBETH COLOR AND LUMINANCE CHART
%
%   The default, scene = sceneCreate, is a Macbeth color checker illuminated
%   by a D65 light source with a mean luminance of 100 cd/m2.  The scene is
%   described only a small number of spatial 64x96 (row,col).  This can be
%   changed using the patchSize argument (default - 16 pixels).  The
%   wavelength  400:10:700 samples, making it efficient to use for experiments.
%
%   Example:
%    scene = sceneCreate('macbeth',32);
%
%    patchSize = 8;
%    spectrum.wave = (380:4:1068)';
%    scene = sceneCreate('macbethEE_IR',patchSize,spectrum);
%
%      {'macbethd65'}  - Create a Macbeth D65 image.  Optional
%         parameter of patch size (default = 16 pixels).
%      {'macbethd50'}         - D50 illuminant
%      {'macbethillc'}        - Illuminant C
%      {'macbethfluorescent'} - Fluorescent illuminant
%      {'macbethtungsten'}    - Tungsten illuminant
%      {'macbethEE_IR'}       - Equal energy extends out to the IR
%
%   You can use sceneAdjustIlluminant() to change the scene SPD.  This
%   function runs on MCC and any other scene with an illuminant.
%
%   The size of the individual patches and the wavelength sampling are both
%   parameters. They can be set using the calling procedure
%
%         patchSizePixels = 16;
%         spectrum.wave = [380:5:720];
%         scene = sceneCreate('macbethTungsten',patchSizePixels,spectrum);
%
%     {L*-steps} - Vertical bars spaced in equal L* steps (dark -> light)
%     scene = sceneCreate('LSteps',barWidth=20,nBars=10,deltaE=10);
%
% REFLECTANCE SAMPLE CHART
%
%   {'reflectance chart'}       - Specify random reflectance samples from
%                                 database. There is always a gray strip at
%                                 the right.   Uses sceneReflectanceChart
%
%       pSize = 24;          % Patch size in pixels
%       sSamples = [64 64];  % Surface samples from the files
%       sFiles{1} = fullfile(isetRootPath,'data','surfaces','reflectances','MunsellSamples_Vhrel.mat');
%       sFiles{2} = fullfile(isetRootPath,'data','surfaces','reflectances','Food_Vhrel.mat');
%       sceneCreate('reflectance chart',pSize,sSamples,sFiles); 
%
% SPATIAL TEST PATTERNS:
%
%      {'ringsRays'}            - Resolution pattern
%      {'harmonic'}             - Harmonics (can be sums of harmonics)
%      {'sweepFrequency'}       - Increasing frequency to the right,
%               increasing contrast upward
%      {'lined65'}              - Line with
%      {'lineee''}              - Line with equal energy spectrum
%      {'pointArray'}           - Point array
%      {'gridlines'}            - Grid lines
%      {'checkerboard'}         - Checkerboard with equal photon spectrum
%      {'frequencyOrientation'} - Demosaicking test pattern, equal photon spectrum
%      {'slantedEdge'}  - Used for ISO spatial resolution, equal photon spectrum
%      {'zonePlate'}   - Circular zone plot, equal photon spectrum
%      {'starPattern'} - Radial lines used to test printers and displays
%
%  Additional parameters are available for several of the patterns.  For
%  example, the harmonic call can set the frequency, contrast, phase,
%  angle, row and col size of the harmonic.  The frequency unit in this
%  case is cycles/image.  To obtain cycles per degree, divide by the field
%  of view.
%
%        parms.freq = 1; parms.contrast = 1; parms.ph = 0;
%        parms.ang= 0; parms.row = 128; parms.col = 128;
%        parms.GaborFlag=0;
%        [scene,parms] = sceneCreate('harmonic',parms);
%
%  See the script s_sceneHarmonics for more examples.  In this example, the
%  illuminant is set so that the mean of the harmonic has a 20%
%  reflectance, like a typical gray card.
%
%  Many of the patterns can have an arbitrary image (row,col) size.  This
%  is possible for whitenoise, impulse1dee,lined65,
%
%         imageSize = 128; scene = sceneCreate('lined65',imageSize);
%
%  Other patterns have different parameters:
%
%         sceneCreate('slantedBar',imageSize,edgeSlope);
%         scene = sceneCreate('checkerboard',pixelsPerCheck,numberOfChecks)
%         scene = sceneCreate('gridlines',imageSize,pixelsBetweenLines);
%         scene = sceneCreate('pointarray',imageSize,pixelsBetweenPoints);
%
% NOISE ANALYSIS TEST PATTERNS
%
%      {'linearIntensityRamp'}  -
%      {'uniformEqualEnergy'}   - Equal energy
%      {'uniformEqualPhoton'}   - Equal photon density
%      {'uniformd65'}           -
%      {'whitenoise'}           - Noise pattern for testing
%
%    The uniform patterns are small by default (32,32).  If you would like
%    them at a higher density (not much point), you can use
%        sceneCreate('uniformD65',256)
%    where 256 is the image size in pixels.
%
% SCENES FROM IMAGES
%      It is possible to create scenes  using data in image fiels.
%      The high dynamic range, multispectral image data included in
%      data\images\MultiSpectral directory are an important source of
%      data.  It is also possible to simply read a tiff or jpeg file and
%      create a scene structure.  These image-based scenes created by the
%      call sceneFromFile which, in turn, calls this function.
%
% Copyright ImagEval Consultants, LLC, 2003.

if ieNotDefined('sceneName'), sceneName = 'default'; end

% Identify the object type
scene.type = 'scene';

sceneName = ieParamFormat(sceneName);

switch sceneName
    case 'default'
        % The user can make a Macbeth with different patch sizes and
        % wavelength sampling, by calling with additional arguments, such
        % as scene = sceneCreate('macbethd65',16,spectrum);
        scene = sceneDefault(scene,'d65');
    case {'macbeth','macbethd65'}
        % sceneCreate('macbethD65',24);
        scene = sceneDefault(scene,'d65',varargin);
    case {'macbethd50'}
        scene = sceneDefault(scene,'d50',varargin);
    case {'macbethc','macbethillc'}
        scene = sceneDefault(scene,'c',varargin);
    case {'macbethfluorescent','macbethfluor'}
        scene = sceneDefault(scene,'fluorescent',varargin);
    case {'macbethtungsten','macbethtung'}
        scene = sceneDefault(scene,'tungsten',varargin);
    case {'macbethee_ir','macbethequalenergyinfrared'}
        % Equal energy illumination into the IR
        % The way to call this would be
        % patchSize = 16;
        % spectrum.wave = 380:4:1068;
        % scene = sceneCreate('macbethEE_IR',patchSize,spectrum)
        scene = sceneDefault(scene,'ir',varargin);
    case {'reflectancechart'}
        % sceneCreate('reflectance chart',pSize,sSamples,sFiles); 
        % There is always a gray strip at the right.
        
        % Patch size in pixels
        pSize = 24;
        % Default surface files
        sFiles{1} = fullfile(isetRootPath,'data','surfaces','reflectances','MunsellSamples_Vhrel.mat');
        sFiles{2} = fullfile(isetRootPath,'data','surfaces','reflectances','Food_Vhrel.mat');
        % Surface samples from the files
        sSamples = [64 64];

        if isempty(varargin)
        else
            pSize = varargin{1};
            if length(varargin) > 1, sSamples = varargin{2}; end
            if length(varargin) > 2, sFiles = varargin{3}; end
        end
        scene = sceneReflectanceChart(sFiles,sSamples,pSize);

    case {'lstar'}
        scene = sceneSet(scene,'name','L-star');
        bWidth = 20; nBars = 10; deltaE = 10;
        if ~isempty(varargin)
            bWidth = varargin{1};
            if length(varargin)>1, nBars  = varargin{2}; end
            if length(varargin)>2, deltaE = varargin{3}; end
        end
        scene = sceneLstarSteps(scene,bWidth,nBars,deltaE);
        
        % Monochrome,RGB and multispectral add only a little.  Mostly created in sceneFromFile
    case {'monochrome','unispectral'}
        % sceneMonochrome is used for images with only one spectral band.
        scene = sceneMonochrome(scene);
    case {'multispectral','hyperspectral'}
        scene = sceneMultispectral(scene);
    case 'rgb'
        if isempty(varargin), scene = sceneRGB;
        else scene = sceneRGB(varargin{1});end
        
    case {'mackay','rayimage','ringsrays'}
        % Also called the Siemens star pattern
        % radF = 24; imSize = 512;
        % vcAddAndSelectObject(sceneCreate('mackay',radF,imSize));
        % sceneWindow();
        radFreq = 8; sz = 256;
        if length(varargin) >= 1, radFreq = varargin{1}; end
        if length(varargin) >= 2, sz = varargin{2}; end
        if length(varargin) >= 3
            wave  = varargin{3};
            scene = sceneSet(scene,'wave',wave);
        end
        scene = sceneMackay(scene,radFreq,sz);
    case {'harmonic','sinusoid'}
        if isempty(varargin),
            [scene,parms] = sceneHarmonic(scene);
        elseif length(varargin) == 1
            parms = varargin{1};
            [scene,parms] = sceneHarmonic(scene,parms);
        elseif length(varargin) == 2
            parms = varargin{1};
            wave = varargin{2};
            [scene,parms] = sceneHarmonic(scene,parms, wave);
        else
            error('Wrong number of parameters! Input params structure and optional wavelengths.')
        end
    case {'sweep','sweepfrequency'}
        % These are always equal photon type.  Could add a third argument
        % for spectral type.
        % sz = 512; maxF = sz/16; sceneCreate('sweepFrequency',sz,maxF);
        sz = 128; maxFreq = sz/16;
        if length(varargin) >= 1, sz = varargin{1}; end
        if length(varargin) >= 2, maxFreq = varargin{2}; end
        scene = sceneSweep(scene,sz,maxFreq);
    case {'ramp','linearintensityramp','rampequalphoton'}
        if isempty(varargin),  sz = 32;
        else                   sz = varargin{1};
        end
        scene = sceneRamp(scene,sz);
    case {'uniform','uniformee','uniformequalenergy'}   %Equal energy
        if isempty(varargin),  sz = 32;
        else                   sz = varargin{1};
        end
        scene = sceneUniform(scene,'ee',sz);
        
    case {'uniformeespecify'}   % Equal energy, specify waveband
        % scene = sceneCreate('uniformEESpecify',sz,wavelength);
        sz = 32; wavelength = 400:10:700;
        if ~isempty(varargin), sz = varargin{1}; end
        if length(varargin) > 1, wavelength = varargin{2}; end
        scene = sceneSet(scene,'wave',wavelength(:));
        scene = sceneUniform(scene,'ee',sz);
    case {'uniformequalphoton','uniformephoton'}        %Equal photon density
        % sceneCreate('uniformEqualPhoton',128);
        if isempty(varargin)
            sz = 32;
            scene = sceneUniform(scene,'ephoton',sz);
        elseif length(varargin) == 1
            sz = varargin{1};
            scene = sceneUniform(scene,'ephoton',sz);
        elseif length(varargin) == 2
            sz = varargin{1};
            wave = varargin{2};
            scene = sceneUniform(scene,'ephoton',sz, wave);
        else
            error('Wrong number of arguments : looking for size (optional), wavelengths (optional)');
        end
    case 'uniformd65'
        % sceneCreate('uniformEqualPhoton',64);
        % We should include an option for wavelength so that we extend into
        % the IR
        if isempty(varargin),  sz = 32;
        else                   sz = varargin{1};
        end
        scene = sceneUniform(scene,'D65',sz);
    case {'uniformbb'}
        % scene = sceneCreate('uniformBB',64,5000,400:700);
        if ~isempty(varargin)
            if length(varargin) >= 1, sz = varargin{1}; end
            if length(varargin) >= 2, cTemp = varargin{2}; end
            if length(varargin) >= 3, scene = sceneSet(scene,'wave',varargin{3}); end
        else
            sz = 32; cTemp = 5000;
        end
        scene = sceneUniform(scene,'BB',sz,cTemp);
    case {'uniformmonochromatic'}
        % scene = sceneCreate('uniform monochromatic',sz,wavelength);
        
        % Create a uniform, monochromatic image.  Used for color-matching
        % analyses.  Set the peak radiance in photons.
        sz = 128; wavelength = 500;
        if length(varargin) >= 1, wavelength = varargin{1}; end
        if length(varargin) >= 2, sz = varargin{2}; end
        
        scene = sceneSet(scene,'wave',wavelength);
        scene = sceneUniform(scene,'ee',sz);
        
    case {'lined65','impulse1dd65'}
        if isempty(varargin), sz = 64;
        else sz = varargin{1};
        end
        scene = sceneLine(scene,'D65',sz);
    case {'lineee','impulse1dee'}
        % scene = sceneCreate('lineee',128,2);
        % scene = sceneCreate('lineee',128,2,380:4:1068);
        sz = 64; offset = 0;
        if length(varargin) >= 1, sz = varargin{1};     end
        if length(varargin) >= 2, offset = varargin{2}; end
        if length(varargin) == 3
            scene = sceneSet(scene,'wave',varargin{3});
        end
        scene = sceneLine(scene,'equalEnergy',sz,offset);
    case {'lineequalphoton','lineep'}
        sz = 64; offset = 0;
        if length(varargin) >= 1, sz = varargin{1};     end
        if length(varargin) >= 2, offset = varargin{2}; end
        scene = sceneLine(scene,'equalPhoton',sz,offset);
    case {'whitenoise','noise'}
        % sceneCreate('noise',[128 128])
        sz = 128; contrast = 20;
        if length(varargin) >= 1, sz = varargin{1}; end
        if length(varargin) >= 2, contrast = varargin{2}; end

        scene = sceneNoise(scene,sz,contrast);
        scene = sceneSet(scene,'name','white noise');

    case {'pointarray','manypoints'}
        % sceneCreate('pointArray',sz,spacing,spectralType);
        sz = 128; spacing = 16; spectralType = 'ep';
        if length(varargin) >= 1, sz           = varargin{1}; end
        if length(varargin) >= 2, spacing      = varargin{2}; end
        if length(varargin) >= 3, spectralType = varargin{3}; end
        scene = scenePointArray(scene,sz,spacing,spectralType);
        
    case {'gridlines','distortiongrid'}
        % sceneCreate('gridlines',sz,spacing,spectralType);
        sz = 128; spacing = 16; spectralType = 'ep';
        if length(varargin) >= 1, sz           = varargin{1}; end
        if length(varargin) >= 2, spacing      = varargin{2}; end
        if length(varargin) >= 3, spectralType = varargin{3}; end
        scene = sceneGridLines(scene,sz,spacing,spectralType);
        
    case {'checkerboard'}
        period = 16; spacing = 8; spectralType = 'ep';
        if length(varargin) >= 1, period       = varargin{1}; end
        if length(varargin) >= 2, spacing      = varargin{2}; end
        if length(varargin) >= 3, spectralType = varargin{3}; end
        scene = sceneCheckerboard(scene,period,spacing,spectralType);
        
    case {'demosaictarget','freqorientpattern','frequencyorientation','freqorient'}
        %   parms.angles = linspace(0,pi/2,5);
        %   parms.freqs =  [1,2,4,8,16];
        %   parms.blockSize = 64;
        %   parms.contrast = .8;
        % scene = sceneCreate('freqorient',parms);
        if isempty(varargin), scene = sceneFOTarget(scene);
        else
            % First argument is parms structure
            scene = sceneFOTarget(scene,varargin{1});
        end
    case {'slantedbar','iso12233','slantededge'}
        % scene = sceneCreate('slantedEdge',sz, slope, fieldOfView, wave);
        % scene = sceneCreate('slantedEdge',128,1.33);  % size, slope
        % scene = sceneCreate('slantedEdge',128,1.33,[], (380:4:1064));  % size, slope, wave
        barSlope = []; fov = []; wave = []; imSize = [];
        if length(varargin) >= 1, imSize = varargin{1}; end
        if length(varargin) >= 2, barSlope = varargin{2};  end
        if length(varargin) >= 3, fov = varargin{3}; end
        if length(varargin) >= 4, wave = varargin{4}; end
        scene = sceneSlantedBar(scene,imSize,barSlope,fov,wave);
        
    case {'zoneplate'}
        scene = sceneZonePlate(scene,384);
    case {'starpattern','radiallines'}
        % Thin radial lines - Useful for testing oriented blur
        %
        % scene = sceneCreate('starPattern');
        % scene = sceneCreate('starPattern',384);
        imSize = 256; spectralType = 'ep'; nLines = 8;
        if length(varargin) >=1, imSize = varargin{1}; end
        if length(varargin) >=2, spectralType = varargin{2}; end
        if length(varargin) >=3, nLines = varargin{3}; end
        scene = sceneRadialLines(scene,imSize,spectralType,nLines);
        
    otherwise
        error('Unknown scene format.');
end

% Initialize scene geometry, spatial sampling
scene = sceneInitGeometry(scene);
scene = sceneInitSpatial(scene);

% Scenes are initialized to a mean luminance of 100 cd/m2.  The illuminant
% is adjusted so that dividing the radiance (in photons) by the illuminant
% (in photons) produces the appropriate peak reflectance (default = 1).
%
% Also, a best guess is made about one known reflectance.
if checkfields(scene,'data','photons') && ~isempty(scene.data.photons)
    
    if isempty(sceneGet(scene,'knownReflectance')) && checkfields(scene,'data','photons')
        
        nWave = sceneGet(scene,'nWave');
        
        % If there is no illuminant yet, set the illuminant to all 1's.
        if isempty(sceneGet(scene,'illuminantPhotons'))
            scene = sceneSet(scene,'illuminantPhotons',ones(1,nWave));
        end
        
        % There is no knownReflectance, so we set the peak radiance to a
        % reflectance of 0.9.
        v = sceneGet(scene,'peakRadianceAndWave');
        wave = sceneGet(scene,'wave');
        idxWave = find(wave == v(2));
        p = sceneGet(scene,'photons',v(2));
        [tmp,ij] = max2(p);
        v = [0.9 ij(1) ij(2) idxWave];
        scene = sceneSet(scene,'knownReflectance',v);
    end
    
    luminance = sceneCalculateLuminance(scene);
    scene = sceneSet(scene,'luminance',luminance);
    
    % This routine also adjusts the illumination level to be consistent
    % with the reflectance and scene photons.
    scene = sceneAdjustLuminance(scene,100);
end

return;

%---------------------------------------------------
function scene = sceneNoise(scene,sz,contrast)
%% Make a spatial white noise stimulus
% contrast is the standard deviation of the N(0,contrast) noise.
% The noise is shifted to a mean of 0.5, and the level is clipped to a
% minimum of 0. 

if ieNotDefined('sz'), sz = [128,128]; end
if ieNotDefined('contrast'), contrast = 0.20; 
elseif contrast > 1, contrast = contrast/100; 
end

scene = initDefaultSpectrum(scene,'hyperspectral');
wave  = sceneGet(scene,'wave');
nWave = sceneGet(scene,'nwave');

% This is an image with reasonable dynamic range (10x).
d   = randn(sz)*contrast + 1; d = max(0,d);
d65 = vcReadSpectra('D65',wave);
illPhotons = Energy2Quanta(wave,d65);

data = zeros(sz(1),sz(2),nWave);
for ii=1:nWave, data(:,:,ii) = d*illPhotons(ii); end

scene = sceneSet(scene,'illuminantPhotons',illPhotons);

% Allocate space for the (compressed) photons
scene = sceneSet(scene,'cphotons',data);

% By setting the fov here, we will not override the value in
% sceneInitSpatial() when this returns
scene = sceneSet(scene,'fov',1);

return;


%----------------------------------
function scene = sceneDefault(scene,illuminantType,args)
%% Default scene is a Macbeth chart with D65 illuminant and patchSize 16
% pixels.

if ieNotDefined('illuminantType'), illuminantType = 'd65'; end
if ieNotDefined('args'), args = []; end

if (isempty(args) || isempty(args{1})), patchSize = 16;
else patchSize = args{1};
end

% Create the scene variable
scene = sceneSet(scene,'type','scene');
if isempty(args) || length(args) < 2 || isempty(args{2})
    scene = initDefaultSpectrum(scene,'hyperspectral');
else    scene = sceneSet(scene,'spectrum',args{2});
end
spectrum = sceneGet(scene,'spectrum');

switch lower(illuminantType)
    case 'd65'
        scene = sceneSet(scene,'name','Macbeth (D65)');
        lightSource = illuminantCreate('D65',[],100,spectrum);
    case 'd50'
        scene = sceneSet(scene,'name','Macbeth (D50)');
        lightSource = illuminantCreate('D50',[],100,spectrum);
    case 'fluorescent'
        scene = sceneSet(scene,'name','Macbeth (Fluorescent)');
        lightSource = illuminantCreate('Fluorescent',[],100,spectrum);
    case 'c'
        scene = sceneSet(scene,'name','Macbeth (Ill C)');
        lightSource = illuminantCreate('illuminantC',[],100,spectrum);
    case 'tungsten'
        scene = sceneSet(scene,'name','Macbeth (Tungsten)');
        lightSource = illuminantCreate('tungsten',[],100,spectrum);
    case 'ir'
        scene = sceneSet(scene,'name','Macbeth (IR)');
        lightSource = illuminantCreate('equalEnergy',[],100,spectrum);
    otherwise
        error('Unknown illuminant type.');
end

% Default distance in meters.
scene = sceneSet(scene,'distance',1.2);

% Scene magnification is always 1.
% Optical images have other magnifications that depend on the optics.
scene = sceneSet(scene,'magnification',1.0);

% The default patch size is 16x16.
spectrum = sceneGet(scene,'spectrum');
macbethChartObject = macbethChartCreate(patchSize,(1:24),spectrum);
scene = sceneCreateMacbeth(macbethChartObject,lightSource,scene);

return;

%--------------------------------------------------
function scene = sceneMonochrome(scene)
%% Default monochrome - should probably be eliminated
scene = sceneSet(scene,'name','monochrome');

% Probably should be eliminated.

% We set a default spectrum, but this is usually over-ridden, see
% sceneFromFile.
scene = initDefaultSpectrum(scene,'monochrome');

return;

%--------------------------------------------------
function scene = sceneMultispectral(scene)
%% Default multispectral structure

scene = sceneSet(scene,'name','multispectral');
scene = initDefaultSpectrum(scene,'multispectral');

return;

%--------------------------------------------------
function scene = sceneRGB(scene)
%% Prepare a scene for RGB data.

if ieNotDefined('scene'), scene.type = 'scene'; end

scene = sceneSet(scene,'name','rgb');
scene = sceneSet(scene,'type','scene');
scene = initDefaultSpectrum(scene,'hyperspectral');

% Set up an illuminant - but it is not nicely scaled.  And we don't have a
% known reflectance.
wave = sceneGet(scene,'wave');
d65 = vcReadSpectra('D65',wave);
illPhotons = Energy2Quanta(wave,d65);
scene = sceneSet(scene,'illuminantPhotons',illPhotons);

return;

%--------------------------------------------------
function scene = sceneMackay(scene,radFreq,sz)
%% Someone (I think Chris Tyler) told me the ring/ray pattern is also called
% the Mackay chart.
%   Reference from Joyce here:
%
% Some people call it the Siemens Star pattern (Wueller).
%
% We fill the central circle with a masking pattern.  The size of the
% central region is at the point when the rays would start to alias.  The
% the circumference of the central circle is 2*pi*r (with r in units of
% pixels).  When the radial frequency is f, we need a minimum of 2f pixels
% on the circumference.  So the circumference is 2*pi*r, so that we want
% the radius to be at least r = f/pi.  In practice that is too exact for
% the digital domain.  So we double the radius.
%

if ieNotDefined('radFreq'), radFreq = 8; end
if ieNotDefined('sz'),      sz = 256; end

scene = sceneSet(scene,'name','mackay');

if ~isfield(scene,'spectrum')
    scene = initDefaultSpectrum(scene,'hyperspectral');
end
nWave = sceneGet(scene,'nwave');

img = imgMackay(radFreq,sz);

% Insert central circle mask
r = round(2*radFreq/pi);  % Find the radius for the central circle

% Find the distance from the center of the image
[X,Y] = meshgrid(1:sz,1:sz); X = X - mean(X(:)); Y = Y - mean(Y(:));
d = sqrt(X.^2 + Y.^2);

% Everything with a distance less than 2r set to mean gray (128) for now.
l = (d < r);
img(l) = 128;  % figure; imagesc(img)

scene = sceneSet(scene,'cphotons',repmat(img,[1,1,nWave]));

% Set up an illuminant
wave = sceneGet(scene,'wave');
illPhotons = ones(size(wave))*sceneGet(scene,'data max');
scene = sceneSet(scene,'illuminantPhotons',illPhotons);

return;

%--------------------------------------------------
function scene = sceneSweep(scene,sz,maxFreq)
%%  These are always equal photon

if ieNotDefined('sz'), sz = 128; end
if ieNotDefined('maxFreq'), maxFreq = sz/16; end

scene = sceneSet(scene,'name','sweep');
scene = initDefaultSpectrum(scene,'hyperspectral');
nWave = sceneGet(scene,'nwave');

img = imgSweep(sz,maxFreq);

scene = sceneSet(scene,'cphotons',repmat(img,[1,1,nWave]));

% Set up an illuminant - but it is not nicely scaled.  And we don't have a
% known reflectance.
wave = sceneGet(scene,'wave');

illPhotons = ones(size(wave))*sceneGet(scene,'data max');
scene = sceneSet(scene,'illuminantPhotons',illPhotons);

return;

%--------------------------------------------------
function [scene,p] = sceneHarmonic(scene,parms, wave)
%% Create a scene of a (windowed) harmonic function.
%
% Harmonic parameters are: parms.freq, parms.row, parms.col, parms.ang
% parms.ph, parms.contrast
%
% Missing default parameters are supplied by imageHarmonic.
%
% The frequency is with respect to the image (cyces/image).  To determine
% cycles/deg, use cpd: freq/sceneGet(scene,'fov');
%

scene = sceneSet(scene,'name','harmonic');

if ieNotDefined('wave')
    scene = initDefaultSpectrum(scene,'hyperspectral');
else
    scene = initDefaultSpectrum(scene, 'custom',wave);
end

nWave = sceneGet(scene,'nwave');

% TODO: Adjust pass the parameters back from the imgHarmonic window. In
% other cases, they are simply attached to the global parameters in
% vcSESSION.  We can get them by a getappdata call in here, but not if we
% close the window as part of imageSetHarmonic
if ieNotDefined('parms')
    global parms;
    h   = imageSetHarmonic; waitfor(h);
    img = imageHarmonic(parms);
    p   = parms;
    clear parms;
else
    [img,p] = imageHarmonic(parms);
end

% To reduce rounding error problems for large dynamic range, we set the
% lowest value to something slightly more than zero.  This is due to the
% ieCompressData scheme.
img(img==0) = 1e-4;
scene = sceneSet(scene,'cphotons',repmat(img,[1,1,nWave]));

% The image mean is always 1.  We treat this as a 20% gray, and therefore
% we set the illuminant to be 5.
illPhotons = 5*ones(1,nWave);
scene = sceneSet(scene,'illuminant photons',illPhotons);

return;

%--------------------------------------------------
function scene = sceneRamp(scene,sz)
%% Intensity ramp (see L-star chart for L* steps)

if ieNotDefined('sz'), sz = 128; end

scene = sceneSet(scene,'name','ramp');
scene = initDefaultSpectrum(scene,'hyperspectral');
nWave = sceneGet(scene,'nwave');

img = imgRamp(sz);

scene = sceneSet(scene,'cphotons',repmat(img,[1,1,nWave]));

return;

%--------------------------------------------------
function scene = sceneUniform(scene,spectralType,sz,varargin)
%% Create a spatially uniform scene.
%
% Various spd types are supported, including d65, blackbody, equal energy,
% equal photon
%

if ieNotDefined('scene'), error('Scene required.'); end
if ieNotDefined('spectralType'), spectralType = 'ep'; end
if ieNotDefined('sz'), sz = 32; end
scene = sceneSet(scene,'name',sprintf('uniform-%s',spectralType));

if isempty(varargin)
    if ~isfield(scene,'spectrum')
        scene = initDefaultSpectrum(scene,'hyperspectral');
    end
end
wave  = sceneGet(scene,'wave');
nWave = sceneGet(scene,'nwave');

d = ones(sz,sz,nWave);

switch lower(spectralType)
    case {'d65'}
        illPhotons = Energy2Quanta(wave,vcReadSpectra('D65',wave));
        scene = sceneSet(scene,'illuminantComment','D65');
    case {'ee','equalenergy','eenergy'}
        illPhotons = Energy2Quanta(wave,ones(nWave,1));
        scene = sceneSet(scene,'illuminantComment','Equal Energy');
    case {'ep','equalphoton','ephoton'}
        illPhotons = ones(nWave,1);
        scene = sceneSet(scene,'illuminantComment','Equal Photon');
    case {'bb','blackbody'}
        if isempty(varargin), cTemp = 5000;
        else                  cTemp = varargin{1};
        end
        illPhotons = blackbody(wave, cTemp, 'photons');
        scene = sceneSet(scene,'illuminantComment',sprintf('Blackbody %.0f Kelvin',cTemp));
    otherwise
        error('Unknown spectral type:%s\n',spectralType);
end

% Set illuminant
scene = sceneSet(scene,'illuminantPhotons',illPhotons);

% Create scene photons
for ii=1:nWave, d(:,:,ii) = d(:,:,ii)*illPhotons(ii); end
scene = sceneSet(scene,'cphotons',d);

return;

%--------------------------------------------------
function scene = sceneLine(scene,spectralType,sz,offset)
%% Create a single line scene.  This is used for computing linespreads and
% OTFs.

if ieNotDefined('spectralType'), spectralType = 'ep'; end
if ieNotDefined('sz'),     sz = 64; end
if ieNotDefined('offset'), offset = 0; end

scene = sceneSet(scene,'name',sprintf('line-%s',spectralType));

if ~isfield(scene,'spectrum')
    scene = initDefaultSpectrum(scene,'hyperspectral');
end
wave    = sceneGet(scene,'wave');
nWave   = sceneGet(scene,'nwave');

% Black is more than zero to prevent HDR problem with ieCompressData
linePos = round(sz/2) + offset;
photons = ones(sz,sz,nWave)*1e-4;
photons(:,linePos,:) = 1;

% Figure out a way to do this using sceneSet.
switch lower(spectralType)
    case {'ep','equalphoton','ephoton'}
        % Equal number of photons at every wavelength
        p = ones(size(wave));
        scene = sceneSet(scene,'illuminantPhotons',p);
    case {'ee','equalenergy','eenergy'}
        % Equal energy at every wavelength.  The large scale factor applied
        % to the number of photons is just to produce a reasonable energy
        % level.
        p = Energy2Quanta(wave,ones(nWave,1));
        % Check that the photons have equal energy
        % e = Quanta2Energy(wave,p')
        scene = sceneSet(scene,'illuminantPhotons',p);
    case 'd65'
        % D65 spectra for the line
        d65 = vcReadSpectra('D65',wave);
        p = Energy2Quanta(wave,d65);
        scene = sceneSet(scene,'illuminantPhotons',p);
    otherwise
        error('Unknown uniform field type.');
end

% Probably a newer faster way to do this.
for ii=1:nWave, photons(:,:,ii) = photons(:,:,ii)*p(ii); end
scene = sceneSet(scene,'cphotons',photons);


return;

%------------------------
function scene = sceneRadialLines(scene,imSize,spectralType,nLines)
%% Create a Siemens Star (radial line) scene.
%
%   scene = sceneRadialLines(scene,imSize,spectralType,nLines)
%
% In this test chart the intensities along lines from the center are
% constant. Measuring on a circle around the center the intensity is a
% harmonic. Hence, frequency varies as a function of radial distance.
%
% Reference:
%   Dieter Wueller thinks this pattern is cool.
%   Digital camera resolution measurement using sinusoidal Siemens stars
%   Proc. SPIE, Vol. 6502, 65020N (2007); doi:10.1117/12.703817
%
% Examples:
%  scene = sceneCreate('radialLines');
%

if ieNotDefined('scene'), error('Scene must be defined'); end
if ieNotDefined('spectralType'), spectralType = 'ep'; end
if ieNotDefined('imSize'), imSize = 256; end
if ieNotDefined('nLines'), nLines = 8; end

scene = sceneSet(scene,'name',sprintf('radialLine-%s',spectralType));
scene = initDefaultSpectrum(scene,'hyperspectral');

% Determine the line angles
radians = pi*(0:(nLines-1))/nLines;
endPoints = zeros(nLines,2);
for ii=1:nLines
    endPoints(ii,:) = round([cos(radians(ii)),sin(radians(ii))]*imSize/2);
end
% plot(endPoints(:,1),endPoints(:,2),'o')

img = zeros(imSize,imSize);

% The routine for drawing lines could be better.
for ii=1:nLines
    x = endPoints(ii,1); y = endPoints(ii,2);
    u = -x; v = -y;
    % Flip so x is the lower one
    if x > 0,
        tmp = [x,y]; x = u; y = v; u = tmp(1); v = tmp(2);
    end
    
    if ~isequal(u,x), slope = (y - v) / (u - x);
        for jj=x:0.2:u,
            kk = round(jj*slope);
            img(round(kk + (imSize/2)) + 1, round(jj + (imSize/2)) + 1) = 1;
        end
    else img(:, (imSize/2) + 1) = 1;
    end
end

img = img(1:imSize,1:imSize);
% To reduce rounding error problems for large dynamic range, we set the
% lowest value to something slightly more than zero.  This is due to the
% ieCompressData scheme.
img(img==0) = 1e-4;
% figure; imagesc(img)

% Create the photon image
wave    = sceneGet(scene,'wave');
nWave   = sceneGet(scene,'nwave');
photons = zeros(imSize,imSize,nWave);

% Figure out a way to do this using sceneSet.
switch lower(spectralType)
    case {'ep','equalphoton','ephoton'}
        % Equal number of photons at every wavelength
        img = img*10^16;
        for ii=1:nWave, photons(:,:,ii) = img; end
    case {'ee','equalenergy','eenergy'}
        % Equal energy at every wavelength.  The large scale factor applied
        % to the number of photons is just to produce a reasonable energy
        % level.
        p = Energy2Quanta(wave,10^16*ones(nWave,1));
        % Check that the photons have equal energy
        % e = Quanta2Energy(wave,p')
        for ii=1:nWave, photons(:,:,ii) = p(ii)*img; end
    case 'd65'
        % D65 spectra for the line
        d65 = vcReadSpectra('D65',wave);
        for ii=1:nWave, photons(:,:,ii) = d65(ii)*img; end
    otherwise
        error('Unknown uniform field type.');
end

scene = sceneSet(scene,'cphotons',photons);

return;

%-----------------------
function scene = sceneFOTarget(scene,parms)
%% Frequency/Orientation target

if ieNotDefined('parms'), parms = []; end

scene = sceneSet(scene,'name','FOTarget');
scene = initDefaultSpectrum(scene,'hyperspectral');
nWave = sceneGet(scene,'nwave');

img = FOTarget('sine',parms);

% Prevent dynamic range problem with ieCompressData
img = ieClip(img,1e-4,1);

% This routine returns an RGB image.  We take the green channel and expand
% it
scene = sceneSet(scene,'cphotons',repmat(img(:,:,2),[1,1,nWave]));

%
wave = sceneGet(scene,'wave');
illPhotons = ones(size(wave))*sceneGet(scene,'data max');
scene = sceneSet(scene,'illuminantPhotons',illPhotons);

return;

%-------------------
function scene = sceneCheckerboard(scene,checkPeriod,nCheckPairs,spectralType)
%% Checkerboard

if ieNotDefined('scene'), error('Scene required'); end
if ieNotDefined('checkPeriod'), checkPeriod = 16; end
if ieNotDefined('nCheckPairs'), nCheckPairs = 8; end
if ieNotDefined('spectralType'), spectralType = 'ep'; end

scene = sceneSet(scene,'name',sprintf('Checker-%s',spectralType));
scene = initDefaultSpectrum(scene,'hyperspectral');
wave  = sceneGet(scene,'wave');
nWave = sceneGet(scene,'nwave');

% The dynamic range of the checkerboard is kept to < 10^4 to prevent
% problems with the rounding error
d = checkerboard(checkPeriod,nCheckPairs);

% Prevent ieCompressData problem.
d = ieClip(d,1e-4,1);

switch lower(spectralType)
    case {'d65'}
        spect = vcReadSpectra('D65',wave);
    case {'ee','equalenergy'}
        spect = Energy2Quanta(wave,ones(nWave,1));
    case {'ep','equalphoton'}
        spect = ones(nWave,1);
    otherwise
        error('Unknown spectral type:%s\n',spectralType);
end

img = zeros(size(d,1),size(d,2),nWave);
for ii=1:nWave, img(:,:,ii) = d*spect(ii); end

% This routine returns an RGB image.  We take the green channel and expand
% it
scene = sceneSet(scene,'cphotons',img);

return

%---------------------------------------------------------------
function scene = sceneSlantedBar(scene,imSize,barSlope,fieldOfView,wave)
%%
%  Slanted bar, 2 deg field of view
%  Slope 2.6 (upper left to lower right)
%  Default size:  384
%
% The scene is set to equal photons across wavelength.

if ieNotDefined('imSize'),      imSize = 384; end
if ieNotDefined('barSlope'),    barSlope = 2.6; end
if ieNotDefined('fieldOfView'), fieldOfView = 2; end
if ieNotDefined('wave'),        wave = 400:10:700; end
scene = sceneSet(scene,'name','slantedBar');
scene = sceneSet(scene,'wave',wave);

sceneW = sceneGet(scene,'wave');
nWave  = sceneGet(scene,'nwave');

% Make the image
imSize = round(imSize/2);
[X,Y] = meshgrid(-imSize:imSize,-imSize:imSize);
img = zeros(size(X));

%  y = barSlope*x defines the line.  We find all the Y values that are
%  above the line
list = (Y > barSlope*X );
img( list ) = 1;

% Prevent dynamic range problem with ieCompressData
img = ieClip(img,1e-4,1);

% Create an equal energy scene.
sEnergy = repmat(img,[1,1,nWave]);

% Re-write this some day.
% Convert that scene energy into photons.  Do this by finding
% the scale factors that convert energy to photons and applying them to
% every wave band
ee         = ones(nWave,1);   % Equal energy vector
e2pFactors = Energy2Quanta(sceneW,ee);  % Energy to photon factor
for ii=1:nWave
    sPhotons(:,:,ii) = sEnergy(:,:,ii)*e2pFactors(ii);
end
% Quanta2Energy(sceneW,e2pFactors)   % All ones, if the world is good.

% Set the scene photons and field of view
scene = sceneSet(scene,'cphotons',sPhotons);
scene = sceneSet(scene,'horizontalfieldofview',fieldOfView);

% We assume target is perfectly reflective (white), so the illuminant is
% the equal energy illuminant; that is, the SPD is all due to the
% illuminant
scene = sceneSet(scene,'illuminantEnergy',ee);

return;

%-----------------------
function scene = sceneZonePlate(scene,imSize,fieldOfView)
%% Circular zone plate image
%

if ieNotDefined('imSize'), imSize = 256; end
if ieNotDefined('fieldOfView'), fieldOfView = 4; end

scene = sceneSet(scene,'name','zonePlate');
scene = initDefaultSpectrum(scene,'hyperspectral');
nWave = sceneGet(scene,'nwave');

img = imgZonePlate(imSize);
% Prevent dynamic range problem with ieCompressData
img = ieClip(img,1e-4,1);

scene = sceneSet(scene,'cphotons',repmat(img,[1,1,nWave]));
scene = sceneSet(scene,'horizontalfieldofview',fieldOfView);

return

%-----------------------
function scene = sceneLstarSteps(scene,barWidth,nBars,deltaE)
%% Scene with vertical bars in equal L* steps
%
% scene = sceneCreate('lstar',50,5,10);
% vcAddAndSelectObject(scene); sceneWindow;

scene = initDefaultSpectrum(scene,'hyperspectral');

% Create the Y values that will define the intensities of the spd.  First,
% equal spaced L* values
L = (0:(nBars-1))*deltaE;
LAB = zeros(nBars,3);
LAB(:,1) = L(:);

% Transform them to Y values
C = makecform('lab2xyz');
XYZ = applyCform(LAB,C);
Y = XYZ(:,2); Y = Y/max(Y(:));
% vcNewGraphWin; plot(Y)

% Create equal photons illuminant
nWave = sceneGet(scene,'nwave');
illPhotons = ones(nWave,1);
scene = sceneSet(scene,'illuminant photons',illPhotons);

% Now, make the photon image
photons = ones(128,barWidth*nBars,nWave);
for ii=1:nBars
    start = barWidth*(ii-1) + 1; stop = barWidth*ii;
    for jj=1:nWave
        photons(:,start:stop,jj) = Y(ii)*illPhotons(jj);
    end
end

% The level is scaled on return to a mean of 100 cd/m2.
scene = sceneSet(scene,'cphotons',photons);

return
