function chromaticityPlot(pts,background,nPix)
% Draw points superimposed on an xy chromaticity diagram in current axes
%
%     chromaticityPlot(pts,[background='gray'],[nPix=256])
%
% The general surrounding background can be gray (by default), white or
% black.
%
%  pts        -  xy values of points on the graph
%  background -  image background color 'gray' (default)
%  nPix       -  Spatial resolution.  
%
% Examples:
% Just the background
%  chromaticityPlot;   % Gray background
%
% A point on backgrounds at different colors and resolutions
%  vcNewGraphWin; pts = [.33,.33];
%  chromaticityPlot(pts,'gray',256);
%  chromaticityPlot(pts,'black');
%  chromaticityPlot(pts,'white',384);
%
% Copyright ImagEval 2011

%% Defaults
if ieNotDefined('pts'), pts = []; end
if ieNotDefined('background'), background = 'gray'; end
if ieNotDefined('nPix'), nPix = 256; end

%% Create a mesh grid of points filled with xy values

% Create nPix (x,y) samples between 0 and 1.
x = linspace(0.001,1,nPix);
y = linspace(0.001,1,nPix);

% Set a default value for Y (cd/m2).  This influences the background
% appearance.
Y_val = 40;

% get xy coordinates and the appropriate set of xyY points
[xx,yy] = meshgrid(x,y);
[nRows,nCols] = size(xx);
xy = horzcat(xx(:),yy(:));

% Set the xy values to xyY with Y = 40 cd/m2.
xyY = horzcat(xy,ones(size(xy,1),1)*Y_val);

%% Color the points outside the XYZ locus
% Points outside the locus are the color of the background

wave = 380:5:700;
spectrumLocus = chromaticity(vcReadSpectra('XYZ',wave));

inPoints = inpolygon(xy(:,1),xy(:,2),spectrumLocus(:,1),spectrumLocus(:,2));
% vcNewGraphWin; imagesc(XW2RGBFormat(color_me,nRows,nCols));
% axis xy; axis equal
nOutside = sum(~inPoints);
w = zeros(1,1,3); w(1,1,:) = 1;
backXYZ = srgb2xyz(w);
switch background
    case 'white'
        backXYZ = backXYZ*Y_val;
    case 'black'
        backXYZ = backXYZ*0;
    case 'gray'
        backXYZ = backXYZ*Y_val/2;
    otherwise
        error('Unknown background %s\n',background);
end

%% Compute XYZ and then sRGB
XYZ = xyy2xyz(xyY);
XYZ(~inPoints,:) = repmat(squeeze(backXYZ(:))',nOutside,1);
XYZ = XW2RGBFormat(XYZ,nRows,nCols);
sRGB = xyz2srgb(XYZ);

%% Now plot the points
imagesc(x,y,sRGB); axis xy; axis equal
if ~isempty(pts)
    hold on
    plot(pts(:,1),pts(:,2),'ko');
    hold off
end

% Tidy up
grid on
set(gca,'xlim',[0 0.8],'ylim',[0 0.85])
xlabel('CIE-x'); ylabel('CIE-y');

return

%
% function p = plotChromaticityDiagram(srgb_img,x,y,background)
% % Render the horseshoe and draw spectrum locus
%
% % load and normalize image of xy diagram
% % tmp = load(filename);
% max_val = max(max(srgb_img(:,:,3)));
% rgb_img = srgb_img/max_val;
%
% % Set figure colors according to background
% switch background
%     case 'black'
%         fg = [0.7 0.7 0.7];
%         bg = [0 0 0];
%     otherwise
%         fg = [0.3 0.3 0.3];
%         bg = [1 1 1];
% end
%
% % Render the image
% vcNewGraphWin;
% p = imagesc(x,y,rgb_img.^0.8);
%
% % Set up the axes and such
% axis image;
% set(gcf,'Color', bg)
% set(gca,'Color', bg,'xcolor', fg, 'ycolor', fg,...
%     'xgrid','on','ygrid','on','xlim',[0, 0.75],'ylim',[0,0.85],...
%     'xcolor',fg, 'ycolor', fg,...
%     'XMinorGrid', 'off', 'YMinorGrid', 'off',...
%     'MinorGridLineStyle', ':','GridLineStyle',':',...
%     'Fontsize',14,'Fontweight','normal',...
%     'ydir', 'normal');
% hold on
%
% % Draw the spectrum locus
% wave = 380:700;
% XYZ  = vcReadSpectra('XYZ',wave);
% spectrumLocus = chromaticity(XYZ);
%
% plot(spectrumLocus(:,1),spectrumLocus(:,2),...
%     'Linewidth', 1','Color', fg);
% hold on
% line([spectrumLocus(1,1),spectrumLocus(end,1)],...
%     [spectrumLocus(1,2),spectrumLocus(end,2)],...
%     'Linewidth', 1, 'Color', fg );
%
% xlabel('x', 'Fontsize', 14, 'color', fg);
% ylabel('y', 'Fontsize', 14, 'color', fg);
%
% return
%
% %% Messing around with a quantified method
%
% % Set image resolution.
% % pts = []; nPix = 128; background = 'white';
% x = linspace(0.001,1,nPix);  % 0 - 1 in nPix steps
% y = linspace(0.001,1,nPix);
%
% % Set a value for Y (cd/m2)
% % Y_val = 40;
%
% % get xy coordinates and the appropriate set of xyY points
% [xx,yy] = meshgrid(x,y);
% [r,c]   = size(xx);
% xyY     = zeros(r,c,3);
% xyY(:,:,1) = xx; xyY(:,:,2)= yy; xyY(:,:,3) = 1;
%
% xyY = RGB2XWFormat(xyY);
% XYZ = XW2RGBFormat(xyY2xyz(xyY),r,c);
% sRGB = xyz2srgb(XYZ);
% vcNewGraphWin; imagesc(sRGB)
%
% %[nRows,nCols] = size(xx);
%
% %nPoints = nRows*nCols;
%
% % xy = zeros(nRows,nCols,2);
% % xy(:,:,1) = xx; xy(:,:,2) = yy;
% % xy = reshape(xy,nPoints,2);


