function rec_to_tif(rec, file_name, rot_flag, contrast_scale)
% REC_TO_TIF converts a tomogram reconstruction volume and saves it as a
% TIFF stack.
%
% Created: 12/26/2015
% =======
%
% Modified: 12/26/2015 "Created"
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
% REC_TO_TIF(rec, file_name) saves the tomographic reconstruction rec to a
% TIFF file with name file_name.tif. Prior to saving, the tomogram is
% rotated so that each element of the 3D TIFF stack corresponds to one x-y
% orthogonal view of the reconstruction volume.
%
% REC_TO_TIF(rec, file_name, rot_flag) does the same as the previous
% usage, but the boolean input rot_flag can be used to specify whether the
% reconstruction volume should be rotated prior to saving. If not rotated,
% each element of the 3D TIFF stack corresponds to one x-z slice of the
% reconstruction volume.
%
% REC_TO_TIF(rec, file_name, rot_flag, contrast_scale) does the same as
% the previous usage, but allows the user to clamp the range of the
% reconstruction volume to an interval of values specified in input
% contrast_scale.
%
% Input:
% =====
% rec            - 2D or 3D tomographic reconstruction volume.
%
% file_name      - Name of the TIFF file to be created. 
%
% rot_flag       - (OPTIONAL=true) Specifies whether the reconstruction
%                  volume should be rotated around its x axis prior to
%                  saving.
%
% contrast_scale - (OPTIONAL) [2,1] vector. Specifies a range of values to
%                  which (rec)'s range should be clamped prior to saving.
%
% Output: None
% ======

if nargin < 3
    rot_flag = true;
end
if nargin < 4
    contrast_scale = [min(rec(:)), max(rec(:))];
end

% Rotate so that each element of the TIFF stack is an x-y orthogonal view
% of the reconstruction volume.
if rot_flag
    rec = rot90_3D(rec, 2, 3);
end

% Clamp, if a contrast_scale was specified, then rescale to [0, 1] for
% imwrite.
rec = imprep(rec, contrast_scale);

% Create the TIFF stack
for i = 1:size(rec, 3)
    imwrite(rec(:, :, i), file_name, 'tif', 'WriteMode', 'append', 'Compression', 'none');
end
end
