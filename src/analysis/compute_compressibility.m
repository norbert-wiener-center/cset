function com_vals = compute_compressibility(u, thresholds)
% COMPUTE_COMPRESSIBILITY computes the compressibility of input array u for
% a set of thresholds specified in input vector (thresholds).
%
% Created: 12/26/2015
% =======
%
% Modified: 12/26/2015 "Created."
% ========  12/31/2015 "Changed name."
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
% com_vals = COMPUTE_COMPRESSIBILITY(u, thresholds) computes the relative
% compressibility of array u for each of the relative thresholds specified
% in (thresholds). For a given threshold t, this returns the proportion of
% elements of u with magnitude less than t*max(abs(u(:))).
%
% Input:
% =====
% u          - Input array to calculate the compressibility of.
%
% thresholds - [T,1] vector of threshold values, specified as a proportion
%              in the range (0,1).
%
% Output:
% ======
% com_vals - [T,1] vector of compressibility values.

T = length(thresholds);
com_vals = zeros(T, 1);

% Flatten, scale array values to [0,1] and make them positive, to simplify
% the relative threshold calculations.
u = abs(u(:)) / max(abs(u(:)));

% Number of elements in u.
N = length(u);

% Calculate compressibility for each threshold.
for i = 1:T
    t = thresholds(i);
    com_vals(i) = sum(u > t) / N;
end
end
