function patches = create_xy_patches(varargin)
% CREATE_XY_PATCHES creates one or more whitened x-y image patches from a
% collection of tomograms. This is used for creating a visual comparison of
% the outputs of multiple reconstruction methods from the same data.
%
% Created: 12/18/2015
% =======
%
% Modified: 12/18/2015 "Created."
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
% patches = CREATE_XY_PATCHES(x1, y1, z1, M, N, recs) creates 2D x-y
% patches from each of the tomographic reconstructions in the cell array
% recs. Each patch is taken from the x-y plane with z coordinate z1, and
% spans indices (y1, x1):(y1 + M, x1 + N).
%
% patches = CREATE_XY_PATCHES(x1, y1, z1, z2, M, N, recs) creates 3D x-y
% patches from each of the tomographic reconstructions in the cell array
% recs. Each patch is a 3D stack of 2D x-y image patches with z coordinates
% from z1:z2. Each 2D patch spans indices (y1, x1):(y1 + M, x1 + N). Useful
% for making GIFs.
%
% Input:
% =====
% x1   - Top-left x coordinate of the patches.
%
% y1   - Top-left y coordinate of the patches.
%
% z1   - Smallest (or only) z coordinate of the patches.
%
% z2   - Largest z coordinate of the patches. May be omitted.
%
% M    - Patch height.
%
% N    - Patch width.
%
% recs - [nt,1] cell array containing the 3D tomographic reconstructions
%        from which the patches are extracted.
%
% Output:
% ======
% patches - [nt,1] cell array containing the image patches.

% Six inputs. Create 2D patches.
if nargin == 6
    x1 = varargin{1};
    y1 = varargin{2};
    z1 = varargin{3};
    M = varargin{4};
    N = varargin{5};
    recs = varargin{6};
    
    % Only a single z-slice.
    zrange = z1;
    
% Seven inputs. Create 3D patches.
elseif nargin == 7
    x1 = varargin{1};
    y1 = varargin{2};
    z1 = varargin{3};
    z2 = varargin{4};
    M = varargin{5};
    N = varargin{6};
    recs = varargin{7};
    
    % More than one z-slice
    zrange = z1:z2;
end
    
% Number of tomograms.
nt = length(recs);

% Initialize the output cell array.
patches = cell(nt, 1);

% Extract patches.
for i = 1:nt
    % Current tomogram.
    rec = recs{i};

    % Rotate the tomogram, to simplify patch creation.
    rec = rot90_3D(rec, 2, 1);

    % Whiten the tomogram, to create a more fair visual comparison
    % between multiple reconstruction methods.
    rec = whiten(rec);

    % Add the patch to output (patches).
    patches{i} = rec(y1:(y1+M-1), x1:(x1+N-1), zrange);
end
end
