% s_skinReflectanceEstimations
%
% The original data is skinMatrix created by sampling the cheek in
% hyperspectral data for 71 faces (see the function sceneReflectanceChart.m)
% this is a  71 x 148 matrix (148 wavelength samples)
%
% Step 1: Original Scene
%   Create a reflectance chart from the skin reflectance data
%   see s_sceneReflectanceCharts
%
% Step 2: Compressed Scene
% Get a set of basis functions and coefficients that describe the spectral
% reflectances - save the coefficients and basis functions for use later
% [imgMean, basis, coef] = hcBasis(double(scene.data.photons));
%   Test: Compare the original reflectance chart with the chart
%   reconstructed using the small basis functions
% The compressed scene will be a row x col x nBases matrix

% Step 3: Sensor output
% Use the Original Scene in a simulation of an
% image sensor.  The output of the simulation is Nvalues from the N
% different color channels in the image sensor
%   Note that the light will be incorporated into the definition of the imaging sensor
% The sensor output will be a row x col x nSensors matrix
%
% Step 4: Estimated Scene - here we find a Nsensors x Nbases matrix to predict the
% coefficients of the skin reflectances
% 
%
% Step 5: We can calculate the rmse between
% Original Scene and Compressed Scene
%   how much information do we lose when we use a linear model
% Compressed Scene and Estimated Scene
%   how much further information loss is there due to the sensor
%
 
%%  Step 1: Original Scene

%   Create a reflectance chart from the skin reflectance data
sFiles = cell(1,1);
sFiles{1} = fullfile(isetRootPath,'data','surfaces','reflectances','HyspexSkinReflectance.mat');
sSamples = 64;

% How many row/col spatial samples in each patch (they are square)
pSize = 10;    % Patch size
load(fullfile(isetRootPath,'data','surfaces','reflectances','HyspexSkinReflectance.mat'));
wave = wavelength;      % Whatever is in the file
grayFlag = 0;  % No gray strip
sampling = 'no replacement';
[scene, ~, reflectance] = sceneReflectanceChart(sFiles,sSamples,pSize,wave,grayFlag,sampling);
%vcNewGraphWin; plot(wave,reflectances)

% Show it on the screen
vcAddAndSelectObject(scene); sceneWindow;

%% Step 2: Compressed Scene

% Compute the svd.  reflectance = U * S * basis'
[basis, S, wgts] = svd(reflectance,'econ');
vcNewGraphWin; plot(wave,basis(:,1:4))

nBasis = 20;
estS = diag(S);
estS((nBasis+1):end) = 0;
estS = diag(estS);
estR = basis*estS*wgts';
figure; plot(estR(:),reflectance(:),'.')

meanR = mean(reflectance,2);
reflectance0 = reflectance - repmat(meanR,[1 size(reflectance,2)]);
[basis0, S0, wgts0] = svd(reflectance0,'econ');
vcNewGraphWin; plot(wave,basis0(:,1:3))
hold on;
plot(wave,meanR,'k-')
hold off


% Change this so we do not take out the mean
% reflectance = sceneGet(scene,'reflectance');
% 
% % How many bases to return.  If < 1, then it refers to variance explained.  
% nBasis = .999;   
% [~, basis, coef] = hcBasis(reflectance,'canonical',nBasis); 
% [imgMean, basis2, coef2] = hcBasis(reflectance,'mean svd',nBasis); 
% 
% wave = sceneGet(scene,'wavelength');
% % plot(wave,imgMean,'k');
% % hold on;
% 
% g = vcNewGraphWin; 
% set(g,'name','Mean included')
% plot(wave,basis(:,1),'r')
% hold on;
% plot(wave,basis(:,2),'g')
% plot(wave,basis(:,3),'b')
% plot(wave,basis(:,4),'c')
% hold off
% 
% f = vcNewGraphWin; 
% set(f,'name','Mean removed');
% plot(wave,imgMean,'k')
% hold on;
% plot(wave,basis2(:,1),'r')
% plot(wave,basis2(:,2),'g')
% plot(wave,basis2(:,3),'b')
% plot(wave,basis2(:,4),'c')
% hold off

% close


% lets see how well it does
%    img2 = imageLinearTransform(coef,basis.basis');
%    [img2,r,c] = RGB2XWFormat(img2);
%    img2 = repmat(imgMean(:),1,r*c) + img2';
%    img2 = XW2RGBFormat(img2',r,c);
%    figure(2); imagesc(sum(img2,3)); colormap(gray); axis image
%% Step 3: Sensor output
% make a function based on s_sensorStackedPixels

%% Step 4: Estimated Scene
