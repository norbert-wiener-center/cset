function setup_pool(parallel_flag, n)
% SETUP_POOL is a convenience wrapper for MATLAB's parpool function, used
% to simplify working with parallelized and non-parallelized runs of the
% CS-ET algorithm. Runs parpool with n workers if parallel_flag is true.
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
% SETUP_POOL(parallel_flag, n) creates a MATLAB parallel pool with n
% workers if parallel_flag is true and no current pool with n workers
% exists.
%
% Input:
% =====
% parallel_flag - Boolean variable. If true, creates a parallel pool. This
%                 setup is used in CS-ET scripts to easily toggle between
%                 parallelized and serial runs of the CS-ET algorithm.
%
% n             - Number of workers in the created parallel pool.
%
% Output: None
% ======

if parallel_flag
    % Check to see if a parallel pool already exists
    P = gcp('nocreate');
    
    % Create a parallel pool if none exists.
    if isempty(P)
        parpool(n);
    % Create a new parallel pool if one exists with a number of workers
    % different from n.
    elseif P.NumWorkers ~= n
        delete(P);
        parpool(n);
    end
end
end
