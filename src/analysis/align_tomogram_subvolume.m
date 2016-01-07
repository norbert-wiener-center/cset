function [subvolume, xcorrs] = align_tomogram_subvolume(tom_large, tom_small)
% ALIGN_TOMOGRAM_SUBVOLUME uses cross-correlation to find a subvolume of
% one tomogram (tom_large) that optimally matches another, smaller tomogram
% (tom_small). Used for comparing reconstructions of the same STEM tilt
% series by two different methods, in this case SIRT and CS-ET, as the SIRT
% reconstructions have a slightly larger first dimension than the CS-ET
% reconstructions.
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
% subvolume = ALIGN_TOMOGRAM_SUBVOLUME(tom_large, tom_small) returns a
% subvolume of tom_large with the same dimensions as tom_small, which
% maximizes the correlation between the whitened subvolume and the whitened
% tom_small. The returned (subvolume) is not whitened.
%
% [subvolume, xcorrs] = ALIGN_TOMOGRAM_SUBVOLUME(tom_large, tom_small) does
% the same, but also returns a vector (xcorrs) of cross-correlations
% between the whitened subvolumes and the whitened tom_small.
%
% Input:
% =====
% tom_large - The larger tomogram, size [M1,N,P] where M1 > M2.
%             ALIGN_TOMOGRAM_SUBVOLUME returns a subvolume of tom_large.
%
% tom_small - The smaller tomogram, size [M2,N,P] where M2 < M1.
%
% Output:
% ======
% subvolume - [M2,N,P] subvolume of tom_large which, when whitened,
%             maximizes the correlation with the whitened tom_small out of
%             all subvolumes of tom_large with the same dimensions.
% 
% xcorrs    - [M1-M2+1,1] vector of cross-correlations between all [M2,N,P]
%             whitened subvolumes of tom_large and the whitened tom_small.

M1 = size(tom_large, 1);
M2 = size(tom_small, 1);

% Number of subvolumes to check.
m = M1 - M2 + 1;

% Correlations betweened whitened tom_small and tom_large subvolumes.
xcorrs = zeros(m, 1);

% Track minimum and maximum correlations, used in specifying a halting
% condition for the search.
cmin = 1; % minimum
cmax = 0; % maximum
cmax_idx = 0; % index of maximum correlation

% Flatten and whiten tom_small.
tom_small = whiten(tom_small(:));

% Compute tom_large subvolume cross-correlations.
for i = 1:m
    % Extract a subvolume.
    subvolume = tom_large(i:(i+M2-1), :, :);
    
    % Flatten and whiten the subvolume.
    subvolume = whiten(subvolume(:));
    
    % Compute correlation between subvolume and tom_small.
    xcorrs(i) = tom_small' * subvolume;
    
    % Update min and max correlations found so far, if need be.
    if xcorrs(i) < cmin
        cmin = xcorrs(i);
    end
    if xcorrs(i) > cmax
        cmax = xcorrs(i);
        cmax_idx = i;
    end
    
    % Halting condition: If correlation has been decreasing for two
    % iterations, the difference between the min and max correlation is
    % greater than cdiff, and the difference between xcorrs(i) and cmin is
    % less than cfrac * (cmax - cmin), halt and return the subvolume with
    % the highest correlation found so far.
    cdiff = 0.1;
    cfrac = 0.95;
    % This halting condition was determined empirically and might not work
    % in some cases. If you suspect you're dealing with such a case, change
    % use_halting_condition to false.
    use_halting_condition = true;
    % Halting check.
    if use_halting_condition && ...
       (cmax - cmin > cdiff) && ...
       (xcorrs(i) - cmin < cfrac * (cmax - cmin)) && ...
       (xcorrs(i) < xcorrs(max(i-1, 1))) && ...
       (xcorrs(max(i-1, 1)) < xcorrs(max(i-2, 1)))
        break
    end
end

% Truncate the xcorrs vector if we halted early.
if i < m
    xcorrs = xcorrs(1:i);
end

% Return the subvolume that yielded the largest xcorr value.
subvolume = tom_large(cmax_idx:(cmax_idx + M2 - 1), :, :);

end
