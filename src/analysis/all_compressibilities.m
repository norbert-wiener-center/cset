function [icomp, tvcomp, wcomp] = all_compressibilities(data, offset, thresholds)
% ALL_COMPRESSIBILITIES computes the compressibility of each 2D slice
% (third coordinate) of input 3D array (data) in the identity, TV, and DB8
% wavelet domains at the thresholds specified in input vector (thresholds).
% The input (offset) is used to calculate the proper compressibility in the
% identity domain, in the case that it has a non-zero background value.
%
% Created: 12/19/2015
% =======
%
% Modified: 12/19/2015 "Created."
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
% [icomp, tvcomp, wcomp] = ALL_COMPRESSIBILITIES(data, offset, thresholds)
% computes the compressibility of each 2D slice (third coordinate) of input
% 3D array (data) in the wavelet, TV, and identity domains at the
% thresholds specified in input vector (thresholds). The input (offset) is
% used to calculate the proper compressibility in the identity domain, in
% the case that it has a non-zero background value.
%
% Input:
% =====
% data       - 3D input array. This function calculates the 2D
%              compressibility of each of the P slices data(:, :, 1), ...,
%              data(:, :, P) in the identity, TV, and DB8 wavelet domains.
%
% offset     - Scalar offset value, should be equal to the background value
%              of the input (data).
%
% thresholds - [T,1] vector of threshold values used in the compressibility
%              calculations. Should be a proportion between 0 and 1.
%
% Output:
% ======
% icomp  - [T,P] array of compressibility values in the identity domain.
%
% tvcomp - [T,P] array of compressibility values in the TV domain.
% 
% wcomp  - [T,P] array of compressibility values in the DB8 wavelet
%          domain.

% Number of threshold values.
T = length(thresholds);

% Number of 2D slices to calculate compressibility of.
P = size(data, 3);

% DB8 wavelet scaling function.
db8_mother = daubcqf(8);

% Initialize output arrays.
icomp = zeros(T, P);
tvcomp = zeros(T, P);
wcomp = zeros(T, P);

% Calculate compressibility in each transform domain.
for i = 1:P
    % Identity domain. Offset data elements by background value (offset).
    icomp(:, i) = compute_compressibility(data(:, :, i) - offset, thresholds);
    % TV domain.
    tvcomp(:, i) = compute_compressibility(tv(data(:, :, i)), thresholds);
    % DB8 wavelet domain.
    wcomp(:, i) = compute_compressibility(mdwt(data(:, :, i), db8_mother), thresholds);
end
end
