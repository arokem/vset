function fullName = ieSaveMultiSpectralImage(fullName,mcCOEF,basis,comment,imgMean,illuminant)
%
%  fullName = ieSaveMultiSpectralImage(fullName,coef,basis,comment,imgMean,illuminant)
%
% Save a Matlab data file containing a multi-spectral image. The image is
% created using routines in the pdcsoft multicapture directory.
% 
% Input arguments
%   mcCOEF  - coefficients (RGB format)
%   basis   - basis functions functions
%   comment
%   imgMean - in some cases we remove the mean before creating the coeffs
%   illuminant structure
%     .wave  are wavelengths in nanometers
%     .data  are illuminant as a function of wavelength in energy units
%
%
% The full path to the data is returned in fullname.
%
% The SPD of the data can be derived from the coefficients and basis
% functions using: 
%
%    spd = rgbLinearTransform(mcCOEF,basis');
%
% See also: mcCreateMultispectralBases, CombineExposureColor
%
%EXAMPLE:
%  ieSaveMultiSpectralImage('c:\user\Matlab\data\Tungsten','MacbethChart-hdrs',mcCOEF,basis,basisLights,illuminant,comment)
%
% Copyright ImagEval Consultants, LLC, 2005.


if ieNotDefined('mcCOEF'),     error('Coefficients required');     end
if ieNotDefined('basis'),      error('Basis function required.');  end
if ieNotDefined('comment'),    comment = sprintf('Date: %s\n',date); end %#ok<NASGU>
if ieNotDefined('illuminant'), error('Illuminant in energy units required'); end

if ieNotDefined('fullName')
    fullName = vcSelectDataFile('stayput','w','mat','Save multispectral data file.');
end
 
% Write out the matlab data file with all of the key information needed.
% Sometimes we save out data approximated usingly on the SVD
% Other times, we use a principal component method and have an image mean
%
if ieNotDefined('imgMean'), 
    save(fullName,'mcCOEF','basis','comment','illuminant');
else
    save(fullName,'mcCOEF','basis','imgMean','comment','illuminant');
end

return;


