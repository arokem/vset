function scene = sceneCreateMacbeth(surface,lightSource,scene)
% Create a hyperspectral scene of the Macbeth chart.
%   
%   scene = sceneCreateMacbeth(surface,lightSource,[scene])
%
% The surface reflectances and light source are specified as arguments.
% The color signal is computed (in photons); these values are then
% attached to the scene structure.
%
% Used by the sceneWindow callbacks to create the Macbeth images.
%
% Copyright ImagEval Consultants, LLC, 2005.
% Programming notes:
%   This code is prior to the ISET development.  
%   As it stands, it requires some fields that are not present elsewhere in
%   the code:
%     lightSource.spectrum and surface.spectrum
%     lightSource.data.photons
%     surface.data
%   This routine should be re-written to conform to more modern structures
%   in ISET.
%
if ~checkfields(lightSource,'spectrum'), error('Bad light source description.');
elseif ~checkfields(surface,'spectrum'), error('Bad surface description.');
elseif ~compareFields(lightSource.spectrum,surface.spectrum)
    error('Mis-match between light source and object spectral fields')
end

% Here we have the reflectance functions at the size of the final image. We
% store the known reflectance.
[r,c,w] = size(surface.data);

% This is particularly unclear code because of the way the Macbeth surfaces
% are stored on disk.  That should change.  I believe what is happening is
% we take the light source in photons and we make it the same size and
% shape as the surface data.  Then we multiply point by point.
photons = lightSource.data.photons(:,ones(1,r*c));
photons = reshape(photons,[size(surface.data,3) c r]);
photons = permute(photons,[3 2 1]);
% [photons/(s sr m^2 nm)]

% We compute the product of the surface reflectance and illuminant photons
% here
scene = sceneSet(scene,'cphotons',surface.data .* photons);

% Store the light source
scene = sceneSet(scene,'illuminantPhotons',lightSource.data.photons);

% Determine a bright reflectance location. This logic is used a couple of
% places and should probably become a function.
pw = sceneGet(scene,'peakRadianceAndWave');
p = sceneGet(scene,'photons',pw(2));
[peak loc] = max2(p);
wave = sceneGet(scene,'wave'); 
idxWave = find(wave == pw(2));
v = peak / lightSource.data.photons(wave == pw(2));
scene = sceneSet(scene,'knownReflectance',[v,loc(1),loc(2),idxWave]);

return;
