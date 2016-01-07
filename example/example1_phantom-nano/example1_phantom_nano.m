%% example1_phantom_nano.m
%
%  Created: 01/02/2016
%  =======
%  
%  Modified: 01/02/2016 "Created."
%  ========
%
%  Author: Matthew Guay
%  ======  mguay@math.umd.edu
%          Applied Mathematics & Statistics, and Scientific Computation
%          Department of Mathematics
%          University of Maryland, College Park
%          Copyright (C) 2016
%
%  Usage:
%  =====
%  Creates CS-ET and WBP reconstructions of the simulated nanoparticle
%  phantom used in (Guay et al., 2016).

%% Setup

% Get the nano phantom dataset, add additional subfolders to the search
% path.
setup('phantom-nano');

% Load the nano phantom;
load('nano.mat');

% If the tilt series doesn't already exist, create it.
if exist('tilt.mat', 'file') ~= 2
    disp('Nano phantom tilt series not found. Creating it now.');
    % Load reconstruction data, to get projection angles.
    load('recdata.mat');
    % Create the tilt series.
    tilt = radon3(nano, recdata.theta);
    % Add shot noise and white noise to the projections.
    rng(1);
    mn = mean(tilt(:));
    tilt = tilt * 5500 / mn;
    tilt = mn / 5500 * poissrnd(tilt);
    tilt = add_gaussian_noise(tilt, 100, 2);
    % Save the tilt series.
    save('data/phantom-nano/tilt.mat', 'tilt');
    % Clean up, even though we'll just load them again in a second.
    clear tilt recdata;
end

% Load the nano phantom projections and additional reconstruction data -
% projection angles, background intensity, and reconstruction dimensions.
[projs, recdata] = get_projs();

% Subsampled projections.
projs2 = projs(:, 1:2:end);
projs3 = projs(:, 1:3:end, :);

% Subsampled projection angles.
theta2 = recdata.theta(1:2:end);
theta3 = recdata.theta(1:3:end);

% Create new recdata structs for the subsampled projections.
recdata2 = recdata;
recdata2.theta = theta2;
recdata3 = recdata;
recdata3.theta = theta3;

%% Optimization parameters

% Split-Bregman parameters.
% 1x (27 tilts)
mu = 0.01;
lambda = mu * 10;
gamma = mu * 30;
kappa = mu * 0;
% 2x (14 tilts)
mu2 = 0.02;
% 3x (9 tilts)
mu3 = 0.03;
% Inner and outer loop counts.
n_inner = 20;
n_outer = 10;

% Conjugate gradient parameters.
cg_tol = 1e-4;
max_cgiter = 12;

% If true, use the GPU for Radon transforms if MATLAB can detect a CUDA
% capable device.
gpu_flag = true;

% Only one z-slice, so no need for parallel processing of multiple slices.
parallel_flag = false;

% Parameter set input to cset.
% 1x, noiseless
run1 = cset_parameters(recdata, mu, lambda, gamma, kappa, n_inner, n_outer, parallel_flag, cg_tol, max_cgiter, gpu_flag);
% 2x, noiseless
run2 = cset_parameters(recdata2, mu2, lambda, gamma, kappa, n_inner, n_outer, parallel_flag, cg_tol, max_cgiter, gpu_flag);
% 3x, noiseless
run3 = cset_parameters(recdata3, mu3, lambda, gamma, kappa, n_inner, n_outer, parallel_flag, cg_tol, max_cgiter, gpu_flag);

%% Reconstructions

% CS-ET
% 1x
tic
rec_cset1 = cset(projs, run1);
toc
% 2x
tic
%rec_cset2 = cset(projs2, run2);
toc
% 3x
tic
rec_cset3 = cset(projs3, run3);
toc

% WBP
% Don't bother to time, it'll be super quick.
% 1x
rec_wbp1 = wbp3(projs, run1);
% 3x
rec_wbp2 = wbp3(projs2, run2);
% 6x
rec_wbp3 = wbp3(projs3, run3);

% Make an image of the 1x and 3x reconstructions, alongside the original
% dataset. If you didn't change any of the split-Bregman or conjugate
% gradient parameters in the previous cell, your imshow output should look
% identical to that PNG.

% Only use the center 256x256 squares of each image, to keep the whole
% thing smaller.
c = 128:383;

% Use to create a small border between pictures in the mosaic.
border_y = ones(256, 2);
border_x = ones(2, 3*256+4);



% The full image mosaic. Compare with example1_phantom_nano.png. 
img = [[nano(c,c) border_y rec_cset1(c,c) border_y rec_wbp1(c,c)]; border_x;
       [nano(c,c) border_y rec_cset3(c,c) border_y rec_wbp3(c,c)]];

% Display img.
imshow(img);

