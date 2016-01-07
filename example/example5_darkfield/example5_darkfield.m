%% example5_darkfield.m
%  
%  Created: 09/21/2015
%  =======
%
%  Modified: 09/21/2015 "Created"
%  ========  12/28/2015 "Changed name, added comments."
%
%  Author: Matthew Guay
%  ======  mguay@math.umd.edu
%          Applied Mathematics & Statistics, and Scientific Computation
%          Department of Mathematics
%          University of Maryland, College Park
%          Copyright (C) 2015
%
%  Usage:
%  =====
%  Creates CS-ET and WBP reconstructions of the experimental dark-field
%  dataset used in (Guay et al., 2016). This dataset consists of a
%  single-axis STEM tilt series taken in the dark-field imaging mode,
%  preprocessed using IMOD (http://bio3d.colorado.edu/imod/).
%
%  *NOTE* This example requires a dataset not included in the cset
%  repository. You can download it from [INSERT LINK HERE].

%% Setup

% *IMPORTANT* If you've downloaded the paper-darkfield dataset,
% substitute your own path in here.
data_loc = '../cset-data/';

% Setup, add the 'paper-darkfield' dataset to the search path.
setup('paper-darkfield', data_loc);

% Get the projection and angle data.
[projs, recdata] = get_projs();

% There's an error in this particular dataset - remove the last projection.
projs = projs(:, 1:(end-1), :)/3e4;

% Subsample for undersampled reconstructions.
projs3 = projs(:, 1:3:end, :);
projs6 = projs(:, 1:6:end, :);
% Subsample projection angles, create new recdata structs for subsampled
% reconstructions.
recdata3 = recdata;
recdata3.theta = recdata.theta(1:3:end);
recdata6 = recdata;
recdata6.theta = recdata.theta(1:6:end);

%% Reconstruction parameter setup

%%% CS-ET parameter setup
% Split-Bregman parameters.
mu = 0.001;
lambda = mu * 7;
gamma = mu * 1;
kappa = mu * 8;
n_inner = 24;
n_outer = 3;

% Conjugate gradient parameters.
cg_tol = 1e-4;
max_cgiter = 12;

% If true, use the GPU for Radon transforms.
gpu_flag = true;

% If true, perform multiple 2D reconstructions in parallel
parallel_flag = false;
% Number of parallel pool workers.
pool_size = 26;
% If parallel_flag is true, this will call MATLAB's parpool function.
setup_pool(parallel_flag, pool_size);

% Create CS-ET parameter structs.
run1 = cset_parameters(recdata, mu, lambda, gamma, kappa, n_inner, n_outer, parallel_flag, cg_tol, max_cgiter, gpu_flag);
run3 = cset_parameters(recdata3, mu, lambda, gamma, kappa, n_inner, n_outer, parallel_flag, cg_tol, max_cgiter, gpu_flag);
run6 = cset_parameters(recdata6, mu, lambda, gamma, kappa, n_inner, n_outer, parallel_flag, cg_tol, max_cgiter, gpu_flag);

%% Single-slice CS-ET reconstruction

% Select a slice basically at random.
idx = 1025-987;
proj = projs(:, :, idx);

% Reconstruct that slice of the darklower tomogram volume.
tic
rec = cset(proj, run1);
toc

% Whiten and scale for contrast.
rc = 3.2 * whiten(rec);

% Display image. Compare with example5_darkfield.png - if you didn't change
% any of the split-Bregman or conjugate gradient parameters in the previous
% cell, your imshow output should look identical to that PNG.
imshow(rc, [-0.03, 0.03])

%% Full reconstructions
% Will probably take a while.

%%% CS-ET
% 1x
tic
reccset1 = cset3(projs, run1);
toc
% 3x
tic
reccset3 = cset3(projs3, run3);
toc
% 6x
tic
reccset6 = cset3(projs6, run6);
toc

%%% WBP
% 1x
tic
recwbp1 = wbp3(projs, run1);
toc
% 3x
tic
recwbp3 = wbp3(projs3, run3);
toc
% 6x
tic
recwbp6 = wbp3(projs6, run6);
toc