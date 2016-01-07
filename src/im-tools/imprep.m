function y = imprep(x, lims)
% IMPREP prepares an array x for use with imwrite by rescaling all inputs
% to the range [0,1]. Array values can be clamped to a range specified in
% the input lims.
%
% Created: 12/19/2015
% =======
%
% Modified: 12/19/2015 "Created"
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
% y = IMPREP(x) rescales the values in x to lie in the interval [0,1].
%
% y = IMPREP(x, lims) first clamps the values of x to the interval
% [lims(1), lims(2)], then rescales the clamped array to the interval
% [0,1].
%
% Input:
% =====
% x    - Input array. Values get rescaled to the interval [0,1].
%
% lims - (OPTIONAL) [2,1] vector of values to clamp x's range to before 
%        rescaling to [0,1]. 
%
% Output:
% ======
% y - Rescaled array, with values in the range [0,1].

% If a second argument is supplied, clamp x's range to [lims(1), lims(2)].
if nargin == 2
    x = min(lims(2), max(lims(1), x));
end

y = (x - min(x(:))) / (max(x(:)) - min(x(:)));
end
