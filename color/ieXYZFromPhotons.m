function XYZ = ieXYZFromPhotons(photons,wave)
%Convert photon spectral power distribution into CIE XYZ
%
%   XYZ = ieXYZFromPhotons(photons,wave)
%
%  This routine converts a spectral power distribution in photons to CIE
%  XYZ values.  The routine converts photons into energy and then calls
%  ieXYZFromEnergy.  See the comments about units in that that routine. 
%
% Copyright ImagEval Consultants, LLC, 2003.

XYZ = ieXYZFromEnergy(Quanta2Energy(wave,photons),wave);

return;

