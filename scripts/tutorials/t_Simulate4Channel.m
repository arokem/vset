%% t_Simulate4Channel.m
%
%  Simulate  a 4-channel sensor and render the data
% 
% The script begins by creating a full spectral scene. Then the scene
% radiance data are passed through the optics to create an irradiance image
% at the sensor surface.  
%
% A 4-channel sensor is created with filters containing Gaussian
% transmittance. The sensor catch is simulated, and the 4-channel data are
% then rendered.
%
% Key functions:  sensorColorFilter
%
% To learn more about a particular function, type "help function" as in
% "help sensorCreate".  There are also online descriptions of the functions
% at http://www.imageval.com/public/Products/ISET/Manual/ISET_Functions.htm
%
% Copyright ImagEval Consultants, LLC, 2010.


%% Introduction
%
% To see a more elaborate script that creates scenes and optical images,
% please review s_SimulateSystem.  Here, we just accept the default scene,
% optical image and focus on creating the 4 channel sensor.
%
s_initISET

%% SCENE
scene = sceneCreate;
wave  = sceneGet(scene,'wave');

% It is often useful to visualize the data in the scene window
vcAddAndSelectObject('scene',scene); sceneWindow;

% When the window appears, you can scale the window size and adjust the
% font size as well (Edit | Change Font Size). There are many other options
% in the pull down menus for cropping, transposing, and measuring scene
% properties.
%
%% OPTICS
%
oi = oiCreate;
oi = oiCompute(scene,oi);
vcAddAndSelectObject('oi',oi); oiWindow;

%% SENSOR

% Create a set of color filters with Gaussian spectral transmittances
cfType = 'gaussian'; 
cPos   = 450:50:600; width = ones(size(cPos))*25;
d.data = sensorColorFilter(cfType,wave, cPos, width); 
d.filterNames = {'b','g','y','r'};
d.comment     = 'Four color simulation for s_Simulate4Channel.m';
d.wavelength  = wave;

% Save them in a tmpFilter file
filterFile = fullfile(pwd,'tmpFilter.mat');
ieSaveColorFilter(d,filterFile)

% How the filters are positioned in the 2x2 block.
filterOrder = [1 2 ; 3 4];
pixel = pixelCreate;

% Initialize the sensor structure.  The spectrum is set to match the
% pixel.spectrum entry.
sensor = sensorCreate('Custom', pixel, filterOrder, filterFile);
sensor = sensorSetSizeToFOV(sensor,sceneGet(scene,'fov')*1.1,scene,oi);

sensor = sensorSet(sensor,'Name','Camera-Simulation');

% Now load the IR filter
fullFileName = fullfile(isetRootPath,'data','sensor','irfilters','irFilter_schott.mat');
irData = ieReadColorFilter(wave,fullFileName); 
sensor = sensorSet(sensor,'irFilter',irData);

% We are now ready to compute the sensor image
sensor = sensorCompute(sensor,oi);

% We can view sensor image in the GUI.  Note that the image that comes up
% shows the color of each pixel in the sensor mosaic. Also, please be aware
% that
%  * The Matlab rendering algorithm often introduces unwanted artifacts
% into the display window.  You can resize the window to eliminate these.
%  * You can also set the display gamma function to brighten the appearance
% in the edit box at the lower left of the window.
vcAddAndSelectObject('sensor',sensor); sensorImageWindow;

% There are a variety of ways to quantify these data in the pulldown menus.
% Also, you can view the individual pixel data either by zooming on the
% image (Edit | Zoom) or by bringing the image viewer tool (Edit | Viewer).
%
% Type 'help iexL2ColorFilter' to find out how to convert data from an
% Excel Spread Sheet to an ISET color filter file or a spectral file/
%
% ISET includes a wide array of options for selecting color filters,
% fill-factors, infrared blocking filters, adjusting pixel properties,
% color filter array patterns, and exposure modes. 
%%  PROCESSOR
%
% The image processing pipeline is managed by the fourth principal ISET
% structure, the virtual camera image (vci).  This structure allows the
% user to set a variety of image processing methods, including demosaicking
% and color balancing

vci = vcimageCreate;
vci = imageSet(vci,'name','No Balance');
vci = imageSet(vci,'scaledisplay',1);
vci = imageSet(vci,'renderGamma',0.6);
vci = imageSet(vci,'demosaicmethod','multichannel');
vci = imageSet(vci,'colorconversionmethod','MCC Optimized');
vci = imageSet(vci,'internalColorSpace','XYZ');

vci = imageSet(vci,'colorBalanceMethod','None');

% The default properties use bilinear demosaicking, no color conversion or
% balancing.  The sensor RGB values are simply set to the display RGB
% values.
vci = vcimageCompute(vci,sensor);

% This shouldn't be here ... but the 4 channel processing isn't right yet
o = imageGet(vci,'result'); o = o/max(o(:)); vci = imageSet(vci,'result',o); 

% As in the other cases, we can bring up a window to view the processed
% data, this time a full RGB image.
vcAddAndSelectObject(vci); vcimageWindow

%% Color conversion

% In the sensor window we used the pulldown under 
% Analyze | Color | Color Conversion Matrix 
% to find an optimal transform matrix for the MCC
m = [...
   -0.0528   -0.2091    1.5196
    0.0050    0.3319   -0.0801
   -0.9505    1.0048   -0.0823
    2.4197   -0.1467   -0.0339
    ];
vci = imageSet(vci,'colorConversionTransform',m);

% We set  the other transforms to the identity, so that the product of all
% the transforms is just the one above.
vci = imageSet(vci,'colorBalanceTransform',eye(3,3));
vci = imageSet(vci,'ics2DisplayTransform',eye(3,3));

% We set the vci to not ask any questions, just use the current matrices.
vci = imageSet(vci,'colorConversionMethod','current matrix');

% Compute and show.
vci = vcimageCompute(vci,sensor);
vcAddAndSelectObject(vci); vcimageWindow

%% PROCESSOR Experiments

% This window offers the opportunity to perform many parameter
% changes and to evaluate certain metric properties of the current system.
% Try the pulldown menu item (Analyze | Create Slanted Bar) and then run
% the pulldown menu (Analyze | ISO12233) to obtain a spatial frequency
% response function for the slanted bar image in the ISO standard.
%

%% END



