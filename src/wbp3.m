function recs = wbp3(projs, run)
% WBP3 performs WBP reconstruction of a 3D tomographic volume from
% parallel-beam STEM tilt series measurements. Makes use of MATLAB's Image
% Processing Toolbox's radon and iradon functions.
%
% Created: 09/27/2015
% =======
%
% Modified: 09/27/2015 "Created."
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
% WBP3(projs, run) performs a reconstruction from measurements 'projs'
% using parameters saved in struct 'run'. See the Input section of cset.m
% for a definition of the fields which run must contain. An additional
% optional field in run is 'parallel_flag', defined in the cset3.m Input
% section.
%
% Input:
% =====
% projs - NxTxP array of Radon transform data. P slices, each slice
%         containing T projections of length N.
% run   - Struct containing information to set up the tomogram volume
%         reconstruction.
%
% Output:
% ======
% recs - MxNxP tomographic reconstruction volume. M is the reconstruction
%        depth (z-dimension), N is the reconstruction width (x-dimension),
%        P is the reconstruction length (y-dimension).

% Check to see if run contains a parallel_flag field. Use its value if so,
% or the default (false) if not.
if isfield(run, 'parallel_flag')
    parallel_flag = run.parallel_flag;
else
    parallel_flag = false;
end

% Pull the relevant parameters from input 'run'.
theta = run.theta;
M = run.M;
N = run.N;
rec_start = N/2 - M/2 + 1;
rec_finish = N/2 + M/2;
rec_interval = rec_start:rec_finish;
P = size(projs, 3);

recs = zeros(M, N, P);
if parallel_flag
    parfor i = 1:P
        proj = squeeze(projs(:, :, i));
        rec = iradon(proj, theta, 'linear', 'Ram-Lak', 1, N);
        recs(:, :, i) = rec(rec_interval, :);
    end
else
    for i = 1:P
        proj = squeeze(projs(:, :, i));
        rec = iradon(proj, theta, 'linear', 'Ram-Lak', 1, N);
        recs(:, :, i) = rec(rec_interval, :);
    end
end
end


