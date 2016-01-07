function patches = create_xz_patches(varargin)
% CREATE_XZ_PATCHES creates one or more whitened x-z image patches from a
% collection of tomograms. This is used for creating a visual comparison of
% the outputs of multiple reconstruction methods from the same data.
%
% Created: 12/19/2015
% =======
%
% Modified: 12/19/2015 "Created."
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
% patches = CREATE_XZ_PATCHES(x1, z1, y1, M, N, recs) creates 2D x-z
% patches from each of the tomographic reconstructions in the cell array
% recs. Each patch is taken from the x-z plane with y coordinate y1, and
% spans indices (z1, x1):(z1 + M, x1 + N).
%
% patches = CREATE_XZ_PATCHES(x1, z1, y1, y2, M, N, recs) creates 3D x-z
% patches from each of the tomographic reconstructions in the cell array
% recs. Each patch is a 3D stack of 2D x-z image patches with y coordinates
% from y1:y2. Each 2D patch spans indices (z1, x1):(z1 + M, x1 + N). Useful
% for making GIFs.
%
% Input:
% =====
% x1   - Top-left x coordinate of the patches.
%
% z1   - Top-left z coordinate of the patches.
%
% y1   - Smallest (or only) y coordinate of the patches.
%
% y2   - Largest y coordinate of the patches. May be omitted.
%
% M    - Patch height.
%
% N    - Patch width.
%
% recs - Cell containing the 3D tomographic reconstructions from which the
%        patches are extracted.
%
% Output:
% ======
% patches - Cell array containing the patches, one cell for each cell in
%           input recs.

% Six inputs. Create 2D patches.
if nargin == 6
    x1 = varargin{1};
    z1 = varargin{2};
    y1 = varargin{3};
    M = varargin{4};
    N = varargin{5};
    recs = varargin{6};
    
    % Only a single y-slice.
    yrange = y1;
    
% Seven inputs. Create 3D patches.
elseif nargin == 7
    x1 = varargin{1};
    z1 = varargin{2};
    y1 = varargin{3};
    y2 = varargin{4};
    M = varargin{5};
    N = varargin{6};
    recs = varargin{7};
    
    % More than one z-slice
    yrange = y1:y2;
end
    
% Number of tomograms.
nt = length(recs);

% Initialize the output cell array.
patches = cell(nt, 1);

% Extract patches.
for i = 1:nt
    % Current tomogram.
    rec = recs{i};

    % Whiten the tomogram, to create a more fair visual comparison
    % between multiple reconstruction methods.
    rec = whiten(rec);

    % Add the patch to output (patches).
    patches{i} = rec(z1:(z1+M-1), x1:(x1+N-1), yrange);
end
end
