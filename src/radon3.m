function projs = radon3(data, theta, num_detectors, gpu_flag, ss_level)
% RADON3 Computes the Radon transform (projections) of a 3D dataset at
% angles specified in input theta using the ASTRA toolbox Radon transform
% on successive 2D slices of the input (data). Supports CPU computation as
% well as CUDA GPU computation.
%
% Created: 12/20/2015
% =======
%
% Modified: 12/20/2015 "Created."
% ========  12/31/2015 "Added a GPU availability check."
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
% projs = RADON3(data, theta) computes the Radon transform of each 2D slice
% of input (data) along the third dimension of (data), at each of the
% angles specified in theta.
%
% projs = RADON3(data, theta, num_detectors) does the same, but specifies
% the number of projection detectors (i.e. the width of the Radon
% projections) in num_detectors.
%
% projs = RADON3(data, theta, num_detectors, gpu_flag) does the same, but
% specifies whether to use ASTRA's CUDA tools to perform the computation by
% passing the boolean input gpu_flag.
%
% projs = RADON3(data, theta, num_detectors, gpu_flag, ss_level) does the
% same, but specifies the level of supersampling to use for ASTRA's
% CUDA-driven Radon transform.
%
% Input:
% =====
% data          - [M,N,P] array containing tomogram volume information.
%
% theta         - [T,1] vector of Radon transform projection angles.
%
% num_detectors - (OPTIONAL=size(data,2)) Number of detectors for each 
%                 projection.
%
% gpu_flag      - (OPTIONAL=true) if true, use the ASTRA CUDA tools to
%                 compute Radon transforms.
%
% ss_level      - (OPTIONAL=2) Level of supersampling for the Radon 
%                 transform. Only does something if gpuflag == true.
%
% Output:
% ======
% tilt_series - [num_detectors,T,P] array. [num_detectors,T] Radon 
%               transform data for each of the P dimension-3 slices of 
%               (data).

if nargin < 3
    num_detectors = size(data,2);
end

if nargin < 4
    gpu_flag = true;
end

if nargin < 5
    ss_level = 2;
end

% Make sure a CUDA-capable GPU is available before trying to use one.
gpu_flag = gpu_flag && (gpuDeviceCount() > 0);

T = length(theta);

% Get the dimensions of the data volume.
[M,N,P] = size(data); 

% Convert projection angles to the radian range required by ASTRA.
thetar = pi/180 * theta;
thetar(thetar < -pi/4) = thetar(thetar < -pi/4) + 2*pi;

% Create ASTRA projection geometry and volume geometry.
proj_geom = astra_create_proj_geom('parallel', 1.0, num_detectors, thetar);
vol_geom = astra_create_vol_geom(M, N);

% A projector structure is required for ASTRA's CPU Radon transforms.
if ~gpu_flag
    projector_id = astra_create_projector('linear', proj_geom, vol_geom);
end

% Create 2D sinogram structure.
sinogram_id = astra_mex_data2d('create','-sino', proj_geom, 0);

% Create 2D image structure
im_id = astra_mex_data2d('create', '-vol', vol_geom, 0);

% Create Radon transform configuration. GPU or CPU based on the value of
% gpuflag.
if gpu_flag
    cfg_rt = astra_struct('FP_CUDA');
    cfg_rt.ProjectionDataId = sinogram_id;
    cfg_rt.VolumeDataId = im_id;
    cfg_rt.DetectorSuperSampling = ss_level;
else
    cfg_rt = astra_struct('FP');
    cfg_rt.ProjectionDataId = sinogram_id;
    cfg_rt.VolumeDataId = im_id;
    cfg_rt.ProjectorId = projector_id;
end

% Create ASTRA Radon transform algorithm
rt_id = astra_mex_algorithm('create', cfg_rt);

% Initialize the final 3D tilt series
projs = zeros(num_detectors,T,P);

% Perform Radon transforms.
for i = 1:P
    % Copy the current data slice to ASTRA's Radon transform algorithm.
    astra_mex_data2d('set', im_id, data(:,:,i));
    
    % Compute the Radon transform
    astra_mex_algorithm('run', rt_id);
    
    % Save the output to tilt_series
    projs(:,:,i) = astra_mex_data2d('get', sinogram_id)';
end

% Clean up after ASTRA.
astra_mex_data2d('clear');
astra_mex_algorithm('clear');
end
