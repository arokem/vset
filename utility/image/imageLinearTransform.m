function imT = imageLinearTransform(im,T)
% Apply a linear transformation to an RGB image 
%
%  imT = imageLinearTransform(im,T)
%
% The image data (im) are in an N x M X W format, (e.g., W=3 and RGB). The
% routine applies a right side multiply to the data.  Specifically, if an
% image point is represented by the row vector, p = [R,G,B] the matrix
% transforms each color point, p, to an output vector pT
%
% This routine works with colorTransformMatrix, which provides access to
% various standard color transformation matrices. 
%
% This routine also works with im in the format (N x M x W) and a T matrix
% size (W x K)
%
% Example:
%   Returns an NxMx3 xyz Image
%     T = colorTransformMatrix('lms2xyz');
%     xyzImage = imageLinearTransform(lmsImage,T);
%
%     T = imageGet(vci,'displayspd');
%     spectralImage = imageLinearTransform(lmsImage,T);
%
% See Also: colorTransformMatrix
%
% Copyright ImagEval Consultants, LLC, 2003.

% Save out the image size information
[r,c,w] = size(im);

if size(T,1) ~= w
    error('image/T data sizes are incorrect. If im is RGB, size(T,1) must be 3.');
end

% We reshape the image data into a r*c x w matrix
%
im = RGB2XWFormat(im);

% Then we multiply and reformat. 
imT = im*T;
imT = XW2RGBFormat(imT,r,c);

return
