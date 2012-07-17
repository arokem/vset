%% t_imageProcessor
%
%   Demonstrate control of the image processing routines using the vci
%   (virtual camera image) methods.
%
% Copyright ImagEval Consultants, LLC, 2007.

%% First, create a simple colorful test scene
scene  = sceneCreate('macbeth'); 
oi     = oiCreate;
sensor = sensorCreate;
scene  = sceneSet(scene,'fov',sensorGet(sensor,'fov'));

% Compute the optical image and sensor from the scene.
oi     = oiCompute(scene,oi);
sensor = sensorCompute(sensor,oi);

%% We now have sensor data, and we are ready to adjust properties and
% experiment with the image processing calls.
vci = vcimageCreate;
vci = imageSet(vci,'name','default');
vci = imageSet(vci,'scaledisplay',1);
vci = imageSet(vci,'renderGamma',0.6);

% First, compute with the default properties.  This uses bilinear
% demosaicing, no color conversion or balancing.  The sensor RGB values are
% simply set to the display RGB values.
vci = vcimageCompute(vci,sensor);
vcAddAndSelectObject(vci);
vcimageWindow

%% Next, compute with the color transform set to MCC Optimized.  The
% optimization, however, takes place in sensor coordinates.
vci = imageSet(vci,'colorconversionmethod','MCC Optimized');
vci = imageSet(vci,'name','MCC');
vci = vcimageCompute(vci,sensor);
vcAddAndSelectObject(vci);
vcimageWindow

%% Next, compute with the internal color space set to XYZ.
vci = imageSet(vci,'internalcs','XYZ');
vci = vcimageCompute(vci,sensor);
vci = imageSet(vci,'name','MCC-XYZ');
vcAddAndSelectObject(vci);
vcimageWindow

%% Next, compute with the color balance set to 'Gray World'.
vci = imageSet(vci,'colorbalancemethod','gray world');
vci = vcimageCompute(vci,sensor);
vci = imageSet(vci,'name','MCC-XYZ-GW');
vcAddAndSelectObject(vci);
vcimageWindow

%% Now, illustrate a different demosaic algorithm
vci = imageSet(vci,'demosaicmethod','Adaptive Laplacian');
vci = vcimageCompute(vci,sensor);
vci = imageSet(vci,'name','MCC-XYZ-GW-AL');
vcAddAndSelectObject(vci);
vcimageWindow

%% End script