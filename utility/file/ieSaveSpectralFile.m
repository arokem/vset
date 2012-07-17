function fullpathname = ieSaveSpectralFile(wavelength,data,comment,fullpathname)
% Save a spectral data ISET data file
%
%   fullpathname = ieSaveSpectralFile(wavelength,data,comment,[fullpathname]);
%
%  This routine specifies the format that ISET spectral data are stored.
%  The data can be read by vcReadSpectra at a later time.
%
%     wavelength:  An N-vector of wavelengths
%     data:        A matrix, NxM, of spectral functions in the columns, 
%     comment:     A string with a comment.
%  [fullpathname]: Optional full path name for output file.
%
% Example:
%    ieSaveSpectralFile(wave,cones,'Stockman Fundamentals from XXX');
% 
%    c = 'foo';
%    wavelength = [400,500,600]';
%    variable = [1,1,1]';
%    fullpathname = ieSaveSpectralFile(wavelength,variable,c)
%    data = vcReadSpectra(fullpathname,[400:50:600])
%
% See Also
%   ieSaveColorFilter, ieReadColorFilter
%
% Copyright ImagEval Consultants, LLC, 2003.

if ieNotDefined('data') || ieNotDefined('wavelength')
    error('data and wavelength must be defined');
end
if ieNotDefined('comment'), comment = 'No comment'; end

if length(wavelength) ~= size(data,1)
    errordlg('The row size of the variable data must match the length of the vector wavelength');
end

if ieNotDefined('fullpathname')
    fullpathname = vcSelectDataFile([isetRootPath,filesep,'data'],'w','mat');
    if isempty(fullpathname)
        disp('User canceled');
        return;
    end
end

save(fullpathname,'wavelength','data','comment');

return;
