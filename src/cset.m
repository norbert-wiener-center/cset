function u = cset(varargin)
% CSET creates a 2D image array reconstructed from STEM tomographic data
% (i.e. Radon transform data).
%
% Created: 09/20/2015
% =======
%
% Modified: 09/20/2015 "Created."
% ========  12/31/2015 "Made cg_tol, max_cgiter, and gpu_flag optional 
%                       arguments."
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
% u = CSET(proj, run) performs a reconstruction using parameters saved in a
% struct 'run'. See Input section below for parameter definitions.
%
% u = CSET(proj, b, theta, M, N, mu, lambda, gamma, kappa, n_inner,
% n_outer) (DEPRECATED) performs a reconstruction using direct entry of all
% parameters. Preferred method of input is using a struct.
%
% u = CSET(proj, b, theta, M, N, mu, lambda, gamma, kappa, n_inner,
% n_outer, cg_tol, max_cgiter, gpu_flag) (DEPRECATED) performs the same
% operation as before, but with additional specifications of optional
% parameters.
%
% u = CSET(proj, b, theta, M, N, mu, lambda, gamma, kappa, n_inner,
% n_outer, cg_tol, max_cgiter, gpu_flag, im_flag) (DEPRECATED) performs the
% same operation as before, but shows an image of the reconstruction u
% after each update iteration if im_flag is true.
%
% u = CSET(proj, b, theta, M, N, mu, lambda, gamma, kappa, n_inner,
% n_outer, cg_tol, max_cgiter, gpu_flag, im_flag, anim_flag) (DEPRECATED)
% peforms the same operation as before, but return an MxNx(n_inner*n_outer)
% array containing each update iteration's output if anim_flag is true.
% Good for animations, but creates very large files.
%
% Input:
% =====
% proj       - NxT array of T 2D Radon transform measurements, each size N.
%
% run        - Struct containing fields 'b', 'theta', 'M', 'N', 'mu',
%              'lambda', 'gamma', 'kappa', 'n_inner', 'n_outer', and
%              optionally, 'cg_tol', 'max_cgiter', 'gpu_flag', 'im_flag',
%              'anim_flag', whose uses are identical to the variables with
%              the same name below.
%
% b          - Reconstruction background value. Estimate this a priori.
%
% theta      - Tx1 vector of Radon projection angles, range of [-90o, 90o].
%
% M          - Reconstruction array z-depth (number of rows).
%
% N          - Reconstruction array x-width (number of columns).
%
% mu         - Split-Bregman parameter - controls u's update step size.
%
% lambda     - Split-Bregman parameter - TV regularization weight.
%
% gamma      - Split-Bregman parameter - L1(I) regularization weight.
%
% kappa      - Split-Bregman parameter - L1(DWT) regularization weight.
%
% n_inner    - Number of inner loop iterations.
%
% n_outer    - Number of outer (Bregman) loop iterations.
%
% cg_tol     - (OPTIONAL=1e-4) Conjugate gradient stopping tolerance.
%
% max_cgiter - (OPTIONAL=12) Maximum CG iterations before termination.
%
% gpu_flag   - (OPTIONAL=true) If true, use CUDA Radon and adjoint Radon 
%              transforms. CSET will only use GPU capabilities if MATLAB
%              detects a CUDA-capable device.
%
% im_flag    - (OPTIONAL=false) If true, display u after each update 
%              iteration.
%
% anim_flag  - (OPTIONAL=false) If true, return an MxNx(n_inner*n_outer) 
%              array containing each iteration of u.
%
% Output:
% ======
% u          - If ~anim_flag, an MxN array containing the reconstruction
%              data at the final iteration of the algorithm. If anim_flag,
%              an MxNx(n_inner*n_outer) array containing the reconstruction
%              data after every iteration of the algorithm.

%% Process inputs.

if nargin == 2
    % Input contains a 'run' struct. Wrap in a try/catch in case the struct
    % has improper fields.
    try
        proj = varargin{1};
        run = varargin{2};
        b = run.b;
        theta = run.theta;
        M = run.M;
        N = run.N;
        mu = run.mu;
        lambda = run.lambda;
        gamma = run.gamma;
        kappa = run.kappa;
        n_inner = run.n_inner;
        n_outer = run.n_outer;
        % Check if cg_tol, max_cgiter, gpu_flag, im_flag and anim_flag
        % fields exist. If not, use default values.
        if isfield(run, 'cg_tol')
            cg_tol = run.cg_tol;
        else
            cg_tol = 1e-4;
        end
        if isfield(run, 'max_cgiter')
            max_cgiter = run.max_cgiter;
        else
            max_cgiter = 12;
        end
        if isfield(run, 'gpu_flag')
            gpu_flag = run.gpu_flag;
        else
            gpu_flag = true;
        end
        if isfield(run, 'im_flag')
            im_flag = run.im_flag;
        else
            im_flag = false;
        end
        if isfield(run, 'anim_flag')
            anim_flag = run.anim_flag;
        else
            anim_flag = false;
        end
        
    catch
        error('Input struct is missing a required field.')
    end
elseif nargin >= 11 && nargin <= 16
    % Direct parameter input.
    proj = varargin{1};
    b = varargin{2};
    theta = varargin{3};
    M = varargin{4};
    N = varargin{5};
    mu = varargin{6};
    lambda = varargin{7};
    gamma = varargin{8};
    kappa = varargin{9};
    n_inner = varargin{10};
    n_outer = varargin{11};
    % Check for optional arguments.
    cg_tol = 1e-4;
    max_cgiter = 12;
    gpu_flag = true;
    im_flag = false;
    anim_flag = false;
    if nargin > 11
        cg_tol = varargin{12};
    end
    if nargin > 12
        max_cgiter = varargin{13};
    end
    if nargin > 13
        gpu_flag = varargin{14};
    end
    if nargin > 14
    	im_flag = varargin{15};
    end
    if nargin > 15
        anim_flag = varargin{16};
    end
else
    % Incorrect number of parameters entered.
    error('Incorrect number of input parameters.')
end

% Only use the GPU if MATLAB can detect a CUDA-capable device.
gpu_flag = gpu_flag && (gpuDeviceCount() > 0);

%% Variable initialization

% Initialize u to 0.
u = zeros(M, N);

% Initialize variables necessary for split Bregman iteration.
x = zeros(M,N); % shrunk Dx(u) 
y = zeros(M,N); % shrunk Dy(u)
z = zeros(M,N); % shrunk I(u)
w = zeros(M,N); % shrunk W(u)
bx = zeros(M,N);
by = zeros(M,N);
bz = zeros(M,N);
bw = zeros(M,N);

if anim_flag
    % Record every frame, rename it to u at the very end.
    u_all = zeros(M, N, n_inner * n_outer);
end

% Convert theta to the radian range required by ASTRA.
thetar = pi / 180 * theta;
thetar(thetar < -pi/4) = thetar(thetar < -pi/4) + 2 * pi;

%% ASTRA setup.
% Supersampling level for GPU Radon transforms. Hard-coded because higher
% levels don't offer much improvement, but >1 prevents some checkerboarding
% artifcats.
supersample_level = 2;

% Create projector, projection geometry, and volume geometry.
proj_geom = astra_create_proj_geom('parallel', 1.0, N, thetar);
vol_geom = astra_create_vol_geom(M, N);
if ~gpu_flag
    % Projector setup required only for CPU implementation.
    projector_id = astra_create_projector('linear', proj_geom, vol_geom);
end

% Create the volume array.
rec_id = astra_mex_data2d('create', '-vol', vol_geom, 0);

% Create the sinogram array.
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);

% Forward and backward projector setup. 
if gpu_flag
    % Create the CUDA Radon transform configuration.
    cfg_fp = astra_struct('FP_CUDA');
    cfg_fp.ProjectionDataId = sinogram_id;
    cfg_fp.VolumeDataId = rec_id;
    cfg_fp.DetectorSuperSampling = supersample_level;
    
    % Create the CUDA adjoint Radon transform configuration.
    cfg_bp = astra_struct('BP_CUDA');
    cfg_bp.ProjectionDataId = sinogram_id;
    cfg_bp.ReconstructionDataId = rec_id;
    cfg_bp.PixelSuperSampling = supersample_level;
else
    % Create the CPU Radon transform configuration.
    cfg_fp = astra_struct('FP');
    cfg_fp.ProjectionDataId = sinogram_id;
    cfg_fp.VolumeDataId = rec_id;
    cfg_fp.ProjectorId = projector_id;
    
    % Create the CPU adjoint Radon transform configuration.
    cfg_bp = astra_struct('BP');
    cfg_bp.ProjectionDataId = sinogram_id;
    cfg_bp.ReconstructionDataId = rec_id;
    cfg_bp.ProjectorId = projector_id;
end

% Create the ASTRA Radon transform algorithm from our configuration.
fp_id = astra_mex_algorithm('create', cfg_fp);

% Create the ASTRA adjoint Radon transform algorithm.
bp_id = astra_mex_algorithm('create', cfg_bp);

%% RWT setup
% Specify the wavelet type (DB8).
wavelet_type = daubcqf(8);

%% Data preprocessing
% Remove background contributions from the Radon data.
background = b * ones(M, N);
background_proj = radon_astra(background);
proj = proj - background_proj;
% Keep an initial copy of this proj for Bregman updates.
proj0 = proj;

% Useful intermediate calculation for the Split-Bregman procedure.
mrtf = mu * adjradon_astra(proj0);

%% Main body

% Keep track of each iteration (used when anim_flag is true).
n_iterations = 0;

% Main update loop.
for outer = 1:n_outer
    % Outer loop - Bregman iterations happen here.
    
    % The full n_inner inner loops are unnecessary for the early stages of
    % the outer loop. Reduce them for better performance (but always keep
    % at least 3).
    actual_n_inner = max(3, ceil(outer / n_outer * n_inner));
    
    for inner = 1:actual_n_inner
        % Inner loop - iteratively solve the L2 minimization subproblem
        % using conjugate gradients, and solve the L1 minimization
        % subproblem using shrinkage operations.
        
        % Increment n_iterations
        n_iterations = n_iterations + 1;
        
        %%% min-L2 subproblem to solve: min ||Ku-rhs||_2 over all possible
        % inputs u, where K is the symmetric linear operator defined by the
        % function @uker below.
        
        % rhs is defined in (Goldstein, Osher 2009) p.12 as rhs^k. Equal to
        % mrtf (mu * adjradon(proj)) plus regularization adjoint terms
        % calculated only for regularization terms with non-zero weights.
        rhs = mrtf;
        if lambda > 0
            rhs = rhs + lambda * (Dxt(x - bx) + Dyt(y - by));
        end
        if gamma > 0
            rhs = rhs + gamma * (z - bz);
        end
        if kappa > 0
            rhs = rhs + kappa * midwt(w - bw, wavelet_type);
        end
        
        % Use MATLAB's built-in CG routine.
        [uc,~] = pcg(@uker, rhs(:), cg_tol, max_cgiter);
        
        % Norm of the difference between the old and new u values.
        update_diff = norm(uc - u(:));
        
        % Replace u with its updated value.
        u = reshape(uc, M, N);
        
        if anim_flag
            % Save this iteration of u with the background added back in.
            u_all(:, :, n_iterations) = u + b;
        end
        
        % Check for lack of convergence (may happen if lambda, gamma, or
        % kappa is way too large).
        if max(isnan(u(:))) == 1
            error('Convergence failed, NaN values detected.')
        end
        
        % Performance consideration - if the differences between two
        % succesive u update values is small enough, break from the inner
        % loop early. 
        if update_diff / (M * N) < 1e-7 && outer < n_outer
            break
        end
        
        % When input im_flag is true, display an image of the
        % reconstruction after each inner loop iteration.
        if im_flag
            if outer == 1 && inner == 1 % First-run setup.
                % Get the current figure display preference, then change
                % parameter 'ImshowBorder' to 'loose' so that figure titles
                % are visible.
                ipt_val = iptgetpref('ImshowBorder');
                iptsetpref('ImshowBorder','loose');
                
                % Create the figure and axes used to display u.
                figure
                a1 = axes;
                
                % Choose magnification such that the (maglevel)% * max(M,N)
                % = 1024
                mag_level = max(100, 100 * 512 / max(M,N));
            end
            
            % Display u. For its title, display the current outer loop
            % iteration, inner loop iteration, and update residual.
            imshow((u + b), [min(u(:) + b(:)),max(u(:) + b(:))], ...
                'InitialMagnification', mag_level, 'Parent', a1)
            title([num2str(outer) ' ' num2str(inner) ' ' num2str(update_diff)])
            drawnow
        end
        
        %%% min-L1 subproblem is solved via the shrinkage routines defined
        % in (Goldstein, Osher 2009).
        
        % x, y variable updates require the x- and y-derivatives of u.
        if lambda > 0
            dxu = Dx(u);
            dyu = Dy(u);
        end
        
        % w variable update requires W(u) (i.e. DWT of u).
        if kappa > 0
            wu = mdwt(u, wavelet_type);
        end
        
        % Shrinkage routines.
        % Isotropic TV shrinkage to update x, y.
        if lambda > 0
            [x, y] = shrinkTV(dxu + bx, dyu + by, 1 / lambda);
        end
        
        % Identity transform shrinkage to update z.
        if gamma > 0
            z = shrinkL1(u + bz, 1 / gamma);
        end
        
        % Wavelet transform shrinkage to update w.
        if kappa > 0
            w = shrinkL1(wu + bw, 1 / kappa);
        end
        
        % Update Bregman parameters.
        if lambda > 0
            bx = bx + dxu - x;
            by = by + dyu - y;
        end
        if gamma > 0
            bz = bz + u - z;
        end
        if kappa > 0
            bw = bw + mdwt(u, wavelet_type) - w;
        end
        
    end % End of inner loop.
    
    % Radon transform of u.
    ru = radon_astra(u);
    
    % "Add back the noise" to the data.
    proj = proj + proj0 - ru;
    
    mrtf = mu * adjradon_astra(proj);
end % End of the outer loop.

%% Cleanup.

% Trim excess u_all frames.
if anim_flag
    u_all(:, :, (n_iterations+1):end) = [];
end

% Clear the mex data allocated by ASTRA.
astra_mex_data2d('clear');
astra_mex_algorithm('clear');

% Reset iptpref 'ImshowBorder'.
if im_flag
    iptsetpref('ImshowBorder', ipt_val);
end

% Finish creating output u.
if anim_flag
    u = u_all;
else
    u = u + b;
end  

%% Shared auxiliary functions.
% Some auxiliary functions are included before the final 'end' defining
% function cset, in order that some variables from cset's workspace may be
% automatically shared with the auxiliary functions. If your version of
% MATLAB supports this, the variables sinogram_id, bp_id, rec_id, fp_id, M,
% N, mu, lambda, and gamma should be colored blue-green.

function im_out = adjradon_astra(rad_in)
    % Wrapper around ASTRA's adjoint Radon transform code.
    astra_mex_data2d('set', sinogram_id, rad_in');
    astra_mex_algorithm('run', bp_id);
    im_out = astra_mex_data2d('get', rec_id);
end

function rad_out = radon_astra(im_in)
    % Wrapper around ASTRA's Radon transform code.
    astra_mex_data2d('set', rec_id, im_in);
    astra_mex_algorithm('run', fp_id);
    rad_out = astra_mex_data2d('get', sinogram_id);
    rad_out = rad_out';
end

function Ku = uker(u_c)
    % Calculates the linear, symmetric operator K applied to input u
    % reshaped to single-column form (hence u_c). Doing this avoids
    % explicitly computing an MNxMN system matrix.
uuk = reshape(u_c,M,N);
% Calculate Radon transform
ruk = radon_astra(uuk);
rtru = mu * adjradon_astra(ruk); %%%%%

Ku = reshape(rtru + lambda * (4*uuk - Dx(uuk) - Dxt(uuk) - Dy(uuk) - Dyt(uuk)) + (gamma + kappa) * uuk, M*N, 1);
end
end % end of function cset

%% Other auxiliary functions
% Functions below this line are also auxiliary functions, but need not
% share any variables from cset's workspace, and so are not included before
% the end of function cset.

function d = Dx(u)
    % Forward x-difference operator.
    [rows,cols] = size(u);
    d = zeros(rows,cols);
    d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
    %d(:,1) = u(:,1)-u(:,cols); % Cyclic boundary condition.
    d(:,1) = 0; % 0 boundary condition.
end

function d = Dxt(u)
    % Forward x-difference adjoint operator (backward difference).
    [rows,cols] = size(u);
    d = zeros(rows,cols);
    d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);
    %d(:,cols) = u(:,cols)-u(:,1); % Cyclic boundary condition.
    d(:,cols) = 0; % 0 boundary condition.
end

function d = Dy(u)
    % Forward y-difference operator.
    [rows,cols] = size(u);
    d = zeros(rows,cols);
    d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
    %d(1,:) = u(1,:)-u(rows,:); % Cyclic boundary condition.
    d(1,:) = 0; % 0 boundary condition.
end

function d = Dyt(u)
    % Forward y-difference adjoint operator (backward difference).
    [rows,cols] = size(u);
    d = zeros(rows,cols);
    d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);
    %d(rows,:) = u(rows,:)-u(1,:); % Cyclic boundary condition.
    d(rows,:) = 0; % 0 boundary condition.
end

function [xs,ys] = shrinkTV(x,y,lambda)
    % Isotropic TV shrinkage routine used by the split-Bregman procedure.
    s = sqrt(x.*conj(x)+y.*conj(y));
    ss = s-lambda;
    ss = ss.*(ss>0);

    s = s+(s<lambda);
    ss = ss./s;

    xs = ss.*x;
    ys = ss.*y;
end

function zs = shrinkL1(z, gamma)
    % L1 shrinkage routine used by the split-Bregman procedure.
    zt = abs(z) - gamma;
    zt = zt.*(zt>0);
    zs = sign(z).*zt;
end
