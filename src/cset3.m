function recs = cset3(projs, params)
% CSET3 performs CS-ET reconstruction of a 3D tomographic volume from
% parallel-beam STEM tilt series measurements.
%
% Created: 09/20/2015
% =======
%
% Modified: 09/20/2015 "Created."
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
% CSET3(projs, params) performs a reconstruction from measurements 'projs'
% using parameters saved in struct 'params'. See the Input section of cset.m
% for a definition of the fields which params must contain. An additional
% optional field in params is 'parallel_flag', defined in the cset3.m Input
% section below.
%
% Input:
% =====
% projs - NxTxP array of Radon transform data. P slices, each slice
%         containing T projections of length N.
% params   - Struct containing fields defined in the Input section of cset.m.
%         Additional optional field parallel_flag should contain a boolean
%         value indicating if the CS-ET algorithm should perform the 2D
%         slice reconstructions in parallel or not.
%
% Output:
% ======
% recs - MxNxP tomographic reconstruction volume. M is the reconstruction
%        depth (z-dimension), N is the reconstruction width (x-dimension),
%        P is the reconstruction length (y-dimension).

% Check to see if params contains a parallel_flag field. Use its value if so,
% or the default (false) if not.
if isfield(params, 'parallel_flag')
    parallel_flag = params.parallel_flag;
else
    parallel_flag = false;
end

% The params's M and N values should be divisible by 2^k for some k>=4, as
% this allows for better wavelet decomposition. Calculate how many wavelet
% decomposition levels are available, issue a warning if it's less than 4.
factor_M = factor(params.M);
n_twos_M = sum(factor_M == 2);
factor_N = factor(params.N);
n_twos_N = sum(factor_N == 2);
n_wavelet_levels = min(n_twos_M, n_twos_N);
if n_wavelet_levels < 4
    warning('params.M and params.N should be divisible by 16 for best wavelet regularization.')
end

% Number of 2D slices to reconstruct.
P = size(projs, 3);
recs = zeros(params.M, params.N, P);

if parallel_flag
    % Reconstruct 2D slices in parallel using parfor.
    parfor i = 1:P
        % Call squeeze to remove singleton 3rd dimension.
        proj = squeeze(projs(:, :, i));
        recs(:, :, i) = cset(proj, params);
    end
else
    % Reconstruct 2D slices sequentially.
    for i = 1:P
        % Call squeeze to remove singleton 3rd dimension.
        proj = squeeze(projs(:, :, i));
        recs(:, :, i) = cset(proj, params);
    end
end

end
