function photons = Energy2Quanta(wavelength,energy)
% Convert energy (watts) to number of photons.
%
%   photons = Energy2Quanta(wavelength,energy)
%
% Energy (e.g., watts) is converted to photons.  Typically energy is a
% matrix of column vectors that measure the energy at the sample
% WAVELENGTHS specified in nanometers, say watts/sr/m2/nm if we are dealing
% with radiance.
%
% The columns are different samples, or sometimes different spatial
% positions. The entries down a column represent the energy as a function
% of wavelength. The entries across a row represent the energy at a single
% wavelength for each of the samples.
%
% SEE ALSO
%   Quanta2Energy().  
%
% CAUTION:  
% Unfortunately, the two routines Quanta2Energy and Energy2Quanta take data
% data in rows vs. cols.  In the fullness of time, this will be changed;
% but it will require some effort because they are inserted in so many
% important places. 
%
% Examples:
%   wave = 400:10:700;  
%   in = [blackbody(wave,5000,'energy'),blackbody(wave,6000,'energy')];
%   p = Energy2Quanta(wave,in);  
%   figure; plot(wave,p);
%
% Notice that in the return, out becomes a row vector, consistent with XW
% format. Also, notice that we have to transpose p to make this work.
% Tragic.
%
%   out = Quanta2Energy(wave,p');      % out is a row vector, XW format
%   figure; plot(wave,in,'ro',wave,out','k-') 
%
% Copyright ImagEval Consultants, LLC, 2003.

if size(wavelength,2) ~= 1 && size(wavelength,1) ~= 1
        error('Energy2Quanta:  Wavelength must be a vector');
else    wavelength = wavelength(:);      % Force to column
end

% Fundamental constants.  
h = vcConstants('h');		% Planck's constant [J sec]
c = vcConstants('c');		% speed of light [m/sec]

if ndims(energy) == 3
    [n,m,w] = size(energy);
    if w ~= length(wavelength)
        error('Energy2Quanta:  energy must have third dimension length equal to numWave');
    end
    energy = reshape(energy,n*m,w)';
    photons = (energy/(h*c)) .* (1e-9 * wavelength(:,ones(1,n*m)));
    photons = reshape(photons',n,m,w);
else
    [n,m] = size(energy);
    if n ~= length(wavelength)
        errordlg('Energy2Quanta:  energy must have row length equal to numWave');
    end
    photons = (energy/(h*c)) .* (1e-9 * wavelength(:,ones(1,m)));
end

return
