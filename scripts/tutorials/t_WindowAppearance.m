% s_WindowAppearance
%
% It is possible to change some of the GUI appearance properties from the
% command line.  Here we illustrate how to make windows invisible, change
% the white point, and change the font size.
%
% Copyright ImagEval Consultants, LLC, 2010
% 
%% Opening and closing windows from the command line

% Let's start by running ISET
ISET

% Close the main window using a command that sets the window parameter
% to off, as in 
ieMainW('visible','off')

% Or bring it back up
ieMainW('visible','on')

% You can do the same with other windows (e.g.,
sceneWindow('visible','off')

sceneWindow('visible','on')

%% Setting spectral power distribution that appears white RGB = (1,1,1)

% Here is a scene with equal photon count across the visible wavelengths
sceneEP = sceneCreate('uniformEqualPhoton');
vcAddAndSelectObject(sceneEP);
sceneWindow;

% By default ISET has an equal photon image as white.  We can change from
% equal photon to equal energy using the command
ieSessionSet('whitePoint','ee');

% Bringin up the scene window, the appearance of the equal photon image
% will change to bluish.
sceneWindow;

% If we now set the scene to be equal energy, it will appear white
sceneEE = sceneCreate('uniformEqualEnergy');
vcAddAndSelectObject(sceneEE);
sceneWindow;

% This affects all renderings in the scene window and in the optics window.
% Here is an example of how the Macbeth color checker is affected by
% changing the white point between 'ee' and 'ep'

% First render it with equal energy set to (1,1,1)
sceneMCC = sceneCreate;
vcAddAndSelectObject(sceneMCC);
ieSessionSet('whitePoint','ee')
sceneWindow

% Now render it with equal photon set to (1,1,1)
ieSessionSet('whitePoint','ep')
sceneWindow

% Finally, it is possible to set D65 to (1,1,1)
ieSessionSet('whitePoint','d65')
sceneWindow

%% Changing the font size in all the windows, you can use this

% This is not fully correct yet.  Close - but needs some work
ieSessionSet('deltaFont',0)

% For the scene window, get the figure handle.
sceneF = ieSessionGet('sceneFigure');
dSize  = 0; ieFontChangeSize(sceneF,dSize);
sceneWindow;

dSize  = 2; ieFontChangeSize(sceneF,dSize);
sceneWindow
ieSessionGet('deltaFont')

dSize  = -2; ieFontChangeSize(sceneF,dSize); 
ieSessionGet('deltaFont')
sceneWindow


dSize  = -1; ieFontChangeSize(sceneF,dSize);

