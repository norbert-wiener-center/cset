%% example3_phantom_complex.m
%  
%  Created: 12/29/2015
%  =======
%
%  Modified: 12/29/2015 "Created."
%  ========
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
%  Creates CS-ET and WBP reconstructions of the complex membrane phantom
%  used in (Guay et al., 2016).

%% Setup

% Get the complex phantom dataset, add additional subfolders to the
% search path.
setup('phantom-complex');

% Load the complex phantom dataset.
load('pcomplex.mat');

% If the tilt series doesn't already exist, create it.
if exist('tilt.mat', 'file') ~= 2
    disp('Complex phantom tilt series not found. Creating it now.');
    % Load reconstruction data, to get projection angles.
    load('recdata.mat');
    % Create the tilt series.
    tilt = radon3(pcomplex, recdata.theta); %#ok<NASGU>
    % Save the tilt series.
    save('data/phantom-complex/tilt.mat', 'tilt');
    % Clean up, even though we'll just load them again in a second.
    clear tilt recdata;
end

% Load the complex phantom projections and additional reconstruction data -
% projection angles, background intensity, and reconstruction dimensions.
[projs, recdata] = get_projs();

% Subsampled projections.
projs3 = projs(:, 1:3:end, :);
projs6 = projs(:, 1:6:end, :);

% Subsampled projection angles.
theta3 = recdata.theta(1:3:end);
theta6 = recdata.theta(1:6:end);

% Create new recdata structs for the subsampled projections.
recdata3 = recdata;
recdata3.theta = theta3;
recdata6 = recdata;
recdata6.theta = theta6;

% Make noisy projections
projsn = add_gaussian_noise(projs, 1200, 1492);
projsn3 = projsn(:, 1:3:end, :);
projsn6 = projsn(:, 1:6:end, :);

%% Optimization parameters

% Split-Bregman parameters.
% 1x, noiseless projections.
mu = 0.005;
lambda = mu * 10;
gamma = mu * 20;
kappa = mu * 0.1;
% 3x, noiseless.
mu3 = 0.006;
lambda3 = mu3 * 4;
gamma3 = mu3 * 10;
kappa3 = mu3 * 0.1;
% 6x, noiseless.
mu6 = 0.01;
lambda6 = mu6 * 2;
gamma6 = mu6 * 10;
kappa6 = mu6 * 1;
% Inner and outer loop counts.
n_inner = 20;
n_outer = 30;

% Conjugate gradient parameters.
cg_tol = 1e-4;
max_cgiter = 12;

% If true, use the GPU for Radon transforms if MATLAB can detect a CUDA
% capable device.
gpu_flag = true;

% If true, run multiple 2D reconstructions in parallel
parallel_flag = false;
% Number of parallel pool workers.
pool_size = 26;
% If parallel_flag is true, this will call MATLAB's parpool function.
setup_pool(parallel_flag, pool_size);

% Parameter set input to cset.
% 1x, noiseless
run1 = cset_parameters(recdata, mu, lambda, gamma, kappa, n_inner, n_outer, parallel_flag, cg_tol, max_cgiter, gpu_flag);
% 3x, noiseless
run3 = cset_parameters(recdata3, mu3, lambda3, gamma3, kappa3, n_inner, n_outer, parallel_flag, cg_tol, max_cgiter, gpu_flag);
% 6x, noiseless
run6 = cset_parameters(recdata6, mu6, lambda6, gamma6, kappa6, n_inner, n_outer, parallel_flag, cg_tol, max_cgiter, gpu_flag);

%% Single slice CS-ET reconstruction.

% Third-coordinate index of the slice to test.
idx = 10;

% Noiseless.
proj = projs(:, :, idx);
tic
rec = cset(proj, run1);
toc
% Whiten and scale for contrast.
rc = 3.2 * whiten(rec);

% Noisy.
projn = projsn(:, :, idx);
tic
recn = cset(projn, run1);
toc
% Whiten and scale for contrast.
rcn = 3.2 * whiten(recn);

% The x-z slice of the original phantom data that we just reconstructed.
pslice = 3.2 * whiten(pcomplex(:, :, idx));

% Use to create an image concatenating the phantom, noiseless, and noisy
% reconstructions with a small white border in between them.
border = 0.05 * ones(2, 256);

% The full concatenated image.
img = [pslice; border; rc; border; rcn];

% Display the image. Compare with example3_phantom_complex.png. If you
% didn't change any of the Split-Bregman or conjugate gradient parameters
% in the previous cell, your imshow output should have the same phantom and
% noiseless reconstructions as the PNG, and a very similar noisy
% reconstruction.
imshow(img, [-0.08, 0.05]);

%% Run the WBP reconstructions

% Noiseless
tic
recw = wbp3(projs, run1);
recw3 = wbp3(projs3, run3);
recw6 = wbp3(projs6, run6);
toc

% Noisy
tic
recwn = wbp3(projsn, run1);
recwn3 = wbp3(projsn3, run3);
recwn6 = wbp3(projsn6, run6);
toc

%% Run the CS-ET reconstructions

% Noiseless
tic
recc = cset3(projs, run1);
toc
tic
recc3 = cset3(projs3, run3);
toc
tic
recc6 = cset3(projs6, run6);
toc

% Noisy
tic
reccn = cset3(projsn, run1);
toc
tic
reccn3 = cset3(projsn3, run3);
toc
tic
reccn6 = cset3(projsn6, run6);
toc