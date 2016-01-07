function snapshot(file_name, tom, idx, min_val, max_val)
% SNAPSHOT saves an image of one 2D x-y slice of a 3D tomogram (tom), or
% the entirety of a 2D input (tom), to a PNG file with name given by
% (file_name).
%
%
% Created: 12/28/2015
% =======
%
% Modified: 12/28/2015 "Created."
% ========
%
% Author: Matthew Guay
% ======  mguay@math.umd.edu
%         Applied Mathematics & Statistics, and Scientific Computation
%         Department of Mathematics
%         University of Maryland, College Park
%         Copyright (C) 2015
%
% Usage:
% =====
% SNAPSHOT(file_name, tom, idx) with ndims(tom)==3, rotates the 3D tomogram
% (tom) about its x axis to obtain an x-y orthogonal view of the volume,
% and saves an image of the x-y slice with z-coordinate (idx) to a PNG file
% with name (file_name). If ndims(tom)==2, SNAPSHOT saves the 2D array
% (tom) to a PNG file with name (file_name).
%
% SNAPSHOT(file_name, tom, idx, min_val, max_val) does the same as the
% previous usage, but first clamps the range of the tomogram values to
% [min_val, max_val].
%
% Input:
% =====
% file_name - Name of the PNG file to which SNAPSHOT saves the indicated
%             tomogram slice.
%
% tom       - 2D or 3D tomographic reconstruction. Could be any other data
%             array of the same dimensions, I guess.
%
% idx       - If input (tom) is a 3D array, this specifies the z-coordinate
%             of the x-y slice to save an image of.
%
% min_val   - (OPTIONAL) minimum value to clamp (tom)'s range to.
%
% max_val   - (OPTIONAL) maximum value to clamp (tom)'s range to.
%
% Output: None
% ======

if nargin < 4
    min_val = min(tom(:));
end

if nargin < 5
    max_val = max(tom(:));
end

% Add '.png' to the end of imgname if it's not there
% already.
if length(file_name) < 5 || ~strcmp('.png', file_name(end-3:end))
    file_name = [file_name '.png'];
end

% If input is 2D instead of 3D, don't bother with the rotation thing, just
% save an image of it.
if ndims(tom) == 3
    tom = rot90_3D(tom, 2, 1);
    img = tom(:, :, idx);
end

% Rescale so that the image is saved properly
img = imprep(img, [min_val, max_val]);
imwrite(img, file_name, 'png');
end
