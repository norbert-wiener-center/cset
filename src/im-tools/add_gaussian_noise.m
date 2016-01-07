function xn = add_gaussian_noise(x, snr, seed)
% ADD_GAUSSIAN_NOISE is an auxiliary function to add Gaussian noise to an
% array, with a specified signal-to-noise ratio.
%
% Created: 12/18/2015
% =======
%
% Modified: 12/18/2015 "Created."
% ========  01/01/2016 "Added the option to control RNG seeding."
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
% xn = ADD_GAUSSIAN_NOISE(x, snr) adds Gaussian noise with variance equal
% to ||x||^2 / snr to each component of the input array (x), , where ||x||
% is the vector Euclidean norm.
%
% xn = ADD_GAUSSIAN_NOISE(x, snr, seed) does the same as the previous
% usage, but seeds MATLAB's RNG with input (seed), so that results may be
% easily reproduced.
%
% Input:
% =====
% x    - Array to which the noise is added.
%
% snr  - Desired signal-to-noise ratio.
%
% seed - (OPTIONAL) RNG seed.
%
% Output:
% ======
% xn - Noisy output array.

% If a third input was provided, use it to seed MATLAB's RNG.
if nargin > 2
    rng(seed);
end

% Make sure snr is > 0.
if snr <= 0
    error('Input snr must be greater than 0.');
end

% RMS of x.
xrms = 1/length(x(:)) * x(:)' * x(:);

% Generate N(0,1) Gaussian noise with the same dimensions as x.
noise = randn(size(x));

% SNR calculation is meaningless for a zero signal.
if xrms == 0
    warning('Input signal is all zeros, desired SNR cannot be achieved.');
    alpha = 1;
else
    % Noise scaling parameter, to give the correct SNR.
    alpha = sqrt(xrms / snr);
end

% Add the noise to x.
xn = x + alpha * noise;
end
