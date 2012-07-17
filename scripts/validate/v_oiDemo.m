%% s_oiDemo
%
% Test optical image functions

% Diffraction limited simulation properties
oi = oiCreate;
plotOI(oi,'otf',[],550);
plotOI(oi,'otf',[],450);

% Human optics
oi = oiCreate('human');
plotOI(oi,'psf',[],420);
plotOI(oi,'psf',[],550);
