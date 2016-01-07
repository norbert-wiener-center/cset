function w = whiten(x)
% WHITEN transforms an input array x to have zero mean and unit variance.
%
% Created: 12/28/2015
% =======
%
% Modified: 12/28/2015 "Created"
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
% w = WHITEN(x) transforms (x) to have zero mean and unit variance and
% outputs the transformed array as (w).
%
% Input:
% =====
% x - An array of arbitrary dimensions.
%
% Output:
% ======
% w - A transformed version of input (x) with zero mean and unit variance.

% Subtract the mean of x from x.
w = x - mean(x(:));

% Rescale w to have unit L2 norm = unit variance.
w = w / norm(w(:));
end
