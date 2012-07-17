function scene = sceneInterpolateW(scene,newWave,preserveLuminance)
%Wavelength interpolation for scene image data
%
%    scene = sceneInterpolateW(scene,[newWave],[preserveLuminance=1])
%
% Interpolate the wavelength dimension of a scene. By default, the
% resampled scene has the same mean luminance as the original scene.
%
% Examples:
%   scene = sceneCreate;
%   scene = sceneInterpolateW(scene,[400:10:700]);
%
% Monochromatic scene
%   scene = sceneInterpolateW(scene,550);
%   vcAddAndSelectObject(scene); sceneWindow;
%
% Do not preserve luminance
%   scene = sceneInterpolateW(scene,[400:2:700],0);
%   vcAddAndSelectObject(scene); sceneWindow;
%
% Copyright ImagEval Consultants, LLC, 2003.

if ieNotDefined('preserveLuminance'), preserveLuminance = 1; end
if ieNotDefined('scene'), scene = vcGetSelectedObject('scene');
elseif ~strcmp(sceneGet(scene,'type'),'scene')
    errordlg('sceneInterpolationW structure not a scene!');
end

handles = ieSessionGet('sceneimagehandle');

% Note the current scene properties
row   = sceneGet(scene,'row');
col   = sceneGet(scene,'col');
% nWave = sceneGet(scene,'nwave');
curWave = sceneGet(scene,'wave');

if ieNotDefined('newWave')
    prompt={'Start (nm)','Stop (nm)','Spacing (nm)'};
    def={num2str(curWave(1)),num2str(curWave(end)),num2str(sceneGet(scene,'binwidth'))};
    dlgTitle='Wavelength resampling';
    lineNo=1;
    val =inputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(val), return; end
    
    low = str2double(val{1}); high = str2double(val{2}); skip = str2double(val{3});
    if high > low,       waveSpectrum.wave = low:skip:high;
    elseif high == low,  waveSpectrum.wave = low;     % User made monochrome, so onlyl 1 sample
    else
        ieInWindowMessage('Bad wavelength ordering:  high < low. Data unchanged.',handles,5);
        return;
    end
else
    waveSpectrum.wave = newWave;
end

% May do different things for different scene types.  At present, the
% Macbeth chart has data (images may have data too) so we can go back to
% the original data set and rebuild the scene with those data.  Other
% images, however, we just interpolate and extrapolate with 0s.
% name = lower(sceneGet(scene,'name'));
% if findstr(name,'macbeth')
%     sceneType = 'macbeth';
% else
%     sceneType = 'unknown';
% end

meanL = sceneGet(scene,'meanluminance');

photons = sceneGet(scene,'photons');
illuminantPhotons = sceneGet(scene,'illuminantPhotons');

% We clear the data to save memory space.
scene = sceneClearData(scene);

% We do this trick to be able to do a 1D interpolation. It is fast ... 2d
% is slow.  The RGB2XW format puts the photons in columns by wavelength.
% The interp1 interpolates across wavelength
photons = RGB2XWFormat(photons)';
newPhotons = interp1(curWave,photons,waveSpectrum.wave)';
newPhotons = XW2RGBFormat(newPhotons,row,col);

newIlluminant = interp1(curWave,illuminantPhotons,waveSpectrum.wave)';

scene = sceneSet(scene,'spectrum',waveSpectrum);
scene = sceneSet(scene,'compressedphotons',newPhotons);
scene = sceneSet(scene,'illuminantPhotons',newIlluminant');

% Store the scene luminance
scene = sceneSet(scene,'luminance',sceneCalculateLuminance(scene));

% For broadband scenes, we generally want to preserve the original mean
% luminance (stored in meanL) despite the resampling. In some cases, such
% as extracting a monochrome scene, we might not want to preserve the mean
% luminance.
if preserveLuminance
    scene = sceneAdjustLuminance(scene,meanL);
end

return;

