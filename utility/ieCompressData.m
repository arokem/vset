function [cData,mn,mx] = ieCompressData(data,bitDepth,mn,mx)
% Compress the image data stored in SCENE and OPTICALIMAGE.
%
%  [cData,mn,mx]  = ieCompressData(data,[bitDepth =16],[mn=min(data(:)],[mx=max(data(:)])
%
%   The data are quantized to uint16  (or uint8) spread over the original
%   range of the data.  The range of the data is recorded and stored as
%   well.  The rounding (compression, uint16) formula is generally
%
%      cData =  uint16 (mxCompress * (data - mn)/(mx - mn));
%      where mxCompress = 2^bitDepth - 1
%
%   The user can either specify the min/max or allow the data min and max
%   to be used.
%
%   If mn > mx, an error is returned.
%   If mn == mx, cData =  uint16 (mxCompress * (data - mn));
%      When using the data to determine mn,mx and mn==mx, this means that
%      the return is all zeros.
%
%   This compression is inverted in the program ieUncompressData.
%
% See sceneGet, sceneSet for examples of the usage.  
%
% Copyright ImagEval Consultants, LLC, 2005.

% TODO:
% Possibly, this program and its inverse should be contained within those
% M-files. 

if ieNotDefined('bitDepth'), bitDepth = 16; end
if ieNotDefined('mn'), mn = min(data(:)); end
if ieNotDefined('mx'), mx = max(data(:)); end

mxCompress = (2^bitDepth) - 1;
% if mx > mn*mxCompress
%     warning('ISET:compressPhotonDR','Photon dynamic range (%f) too high for compressed photons \n',mx/mn);
% end

if mn > mx,      error('Min/Max error.');
elseif mn == mx, s = mx;        % This case is potential trouble.
else             s = (mx - mn);
end

% We will handle large arrays a little differently, wavelength by
% wavelength
[r,c,w] = size(data);

switch bitDepth
    case 16
        try
            % All at once for most images
            if w > 31
                cData = zeros(r,c,w,'uint16');
                h = waitbar(0,'Compressing photons');
                for ii=1:w
                    waitbar(ii/w,h);
                    cData(:,:,ii) = uint16(round(mxCompress*(data(:,:,ii) - mn)/(s)));
                end
                close(h)
            else
                cData = uint16(round(mxCompress * (data - mn)/(s)));
            end

        catch
            % Probably not needed any more because new condition added
            % above.  Stay worried for a while.
            cData = zeros(r,c,w,'uint16');
            h = waitbar(0,'Compressing photons');
            for ii=1:w
                waitbar(ii/w,h);
                cData(:,:,ii) = uint16(round(mxCompress*(data(:,:,ii) - mn)/(s)));
            end
            close(h)
        end
        
    case 8
        % Unused, really.
        cData = uint8(round(mxCompress * (data - mn)/(s)));
    otherwise
        error('Unknown bit depth.');
end

return;