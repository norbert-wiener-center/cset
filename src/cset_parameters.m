function params = cset_parameters(varargin)
% CSET_PARAMETERS creates a struct containing all of the parameters needed
% to perform a 2D or 3D CS-ET reconstruction.
%
% Created: 09/21/2015
% =======
%
% Modified: 09/21/2015 "Created."
% ========  12/31/2015 "Added a second input method using the recdata
%                       struct."
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
% params = CSET_PARAMETERS(b, theta, M, N, mu, lambda, gamma, kappa,
% n_inner, n_outer, parallel_flag=, cg_tol=, max_cgiter=, gpu_flag=,
% im_flag=, anim_flag=) creates struct (params) with the first ten
% mandatory fields (b, ..., n_outer), as well as any of the optional flag
% fields provided. Field names are the same as the input names listed here.
%
% params = CSET_PARAMETERS(recdata, mu, lambda, gamma, kappa, n_inner,
% n_outer, parallel_flag=, cg_tol=, max_cgiter=, gpu_flag=, im_flag=,
% anim_flag=) creates the same output struct, but uses the input struct
% (recdata), which already has fields 'b', 'theta', 'M', and 'N', to create
% those same fields in (params). All other inputs are identical to the
% previous usage.
%
% Input:
% =====
% b             - Scalar reconstruction background value.
% theta         - Tx1 vector of measurement angles.
% M             - Reconstruction depth (z-dimension, 1st index).
% N             - Reconstruction width (x-dimension, 2nd index).
% mu            - Split-Bregman step-size parameter.
% lambda        - Split-Bregman TV regularization parameter.
% gamma         - Split-Bregman L1(I) regularization parameter.
% kappa         - Split-Bregman L1(DWT) regularization parameter
% n_inner       - Split-Bregman number of inner loop iterations.
% n_outer       - Split-Bregman number of outer loop iterations.
% parallel_flag - (OPTIONAL=false) If true, run multiple 2D reconstructions
%                 in parallel.
% cg_tol        - (OPTIONAL=1e-4) Conjugate gradient stopping tolerance.
% max_cgiter    - (OPTIONAL=12) Conjugate gradient max iterations before
%                 termination.
% gpu_flag      - (OPTIONAL=true) If true, use CUDA Radon and adjoint Radon
%                 transforms.
% im_flag       - (OPTIONAL=false) If true, display CS-ET reconstruction
%                 after each update iteration.
% anim_flag     - (OPTIONAL=false) If true, cset() returns an array
%                 containing each iteration of the reconstruction. Good for
%                 animations.
%
% Output:
% ======
% params           - Struct containing a field for each of the inputs, each
%                 with the same name and containing the same data.

% Initialize the params struct.
params = struct;

% If the first input is a struct, assume it's the recdata struct with
% fields b, theta, M, and N.
if isstruct(varargin{1})
    data = varargin{1};
    % Wrap in try/catch to throw an error if the data struct doesn't have
    % the correct fields.
    try
        params.b = data.b;
        params.theta = data.theta;
        params.M = data.M;
        params.N = data.N;
    catch
        error('Input struct does not have the correct fields.');
    end
    params.mu = varargin{2};
    params.lambda = varargin{3};
    params.gamma = varargin{4};
    params.kappa = varargin{5};
    params.n_inner = varargin{6};
    params.n_outer = varargin{7};
    if nargin > 7
        params.parallel_flag = varargin{8};
    end
    if nargin > 8
        params.cg_tol = varargin{9};
    end
    if nargin > 9
        params.max_cgiter = varargin{10};
    end
    if nargin > 10
        params.gpu_flag = varargin{11};
    end
    if nargin > 11
        params.im_flag = varargin{12};
    end
    if nargin > 12
        params.anim_flag = varargin{13};
    end

% Otherwise, assume the alternate form of input.
else
    params.b = varargin{1};
    params.theta = varargin{2};
    params.M = varargin{3};
    params.N = varargin{4};
    params.mu = varargin{5};
    params.lambda = varargin{6};
    params.gamma = varargin{7};
    params.kappa = varargin{8};
    params.n_inner = varargin{9};
    params.n_outer = varargin{10};
    if nargin > 10
        params.parallel_flag = varargin{11};
    end
    if nargin > 11
        params.cg_tol = varargin{12};
    end
    if nargin > 12
        params.max_cgiter = varargin{13};
    end
    if nargin > 13
        params.gpu_flag = varargin{14};
    end
    if nargin > 14
        params.im_flag = varargin{15};
    end
    if nargin > 15
        params.anim_flag = varargin{16};
    end
end
end
