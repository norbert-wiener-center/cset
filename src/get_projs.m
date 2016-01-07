function [projs, recdata] = get_projs(projs_name, recdata_name)
% GET_PROJS reads in TIF, MRC, or MAT file named projs_name, along with
% additional required reconstruction data from a mat file named
% recdata_name.
%
% Created: 09/21/2015
% =======
%
% Modified: 09/21/205 "Created."
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
% [projs, recdata] = GET_PROJS(projs_name, recdata_name) reads in TIF, MRC
% or MAT file named projs_name, along with additional data from mat file
% recdata_name. Output projs is a 3D array of NxTxP projection data (P
% projections, each with T measurements each of length N). Output recdata
% is a struct containing five fields: theta, M, N, P, b.
%
% [projs, recdata] = GET_PROJS() reads in TIF, MRC or MAT files with the
% default names: 'tilt.mrc'/'tilt.tif'/'tilt.mat' and 'recdata.mat'. This
% should work if you've already called csetup with one of the official
% datasets.
%
% Input:
% =====
% projs_name   - (OPTIONAL, default='tilt.***') Name of the TIFF or MRC
%                file to read in. If none is provided, the default is
%                called automatically based on which type of file is
%                detected.
% recdata_name - (OPTIONAL, default='recdata.mat') Name of the MAT file
%                containing the recdata struct, which holds other data
%                required for performing a reconstruction. See output
%                recdata's definition below for an explanation of the
%                fields.
%
% Output:
% ======
% projs   - NxTxP array of projection data (P projections, each with T
%           measurements each of length N).
% recdata - Struct containing five fields:
%               .theta - Tx1 vector of measurement angles.
%               .M     - Reconstruction depth (z-dimension).
%               .N     - Reconstruction width (x-dimension).
%               .P     - Reconstruction length (y-dimension).
%               .b     - Reconstruction background value.

% Tilt file type. Should be 'tif', 'mrc', or 'mat' to load properly.
tilt_type = '';

if nargin < 1
    % Use the default projs_name input. Check the search path to see if
    % it's a TIFF, MRC, or MAT file.
    if exist('tilt.tif', 'file') == 2
        projs_name = 'tilt.tif';
        tilt_type = 'tif';
        
    elseif exist('tilt.tiff', 'file') == 2
        projs_name = 'tilt.tiff';
        tilt_type = 'tif';
        
    elseif exist('tilt.mrc', 'file') == 2
        projs_name = 'tilt.mrc';
        tilt_type = 'mrc';
        
    elseif exist('tilt.mat', 'file') == 2
        projs_name = 'tilt.mat';
        tilt_type = 'mat';
        
    else
        % Tilt file was not found.
        error('Tilt series file was not found.');
    end
end

if nargin < 2
    % Use the default recdata_name input.
    recdata_name = 'recdata.mat';
end

% If the user supplied a tilt series file name, check its file type.
if nargin > 1
    if length(projs_name) < 5
        % Any name needs at least 4 characters, including file extension.
        error('Invalid tilt series file name. Did you include the file extension?');
        
    elseif strcmp(projs_name((end-3):end), 'tiff') || strcmp(projs_name((end-2):end), 'tif')
        % Check to see if projs_name is a TIFF file.
        tilt_type = 'tif';
        
    elseif strcmp(projs_name((end-2):end), 'mrc')
        % Check to see if projs_name is an MRC file.
        tilt_type = 'mrc';
        
    elseif strcmp(projs_name((end-2):end), 'mat')
        % Check to see if projs_name is a MAT file.
        tilt_type = 'mat';
        
    else
        error('Invalid tilt series file name extension.');
    end
end

% Read in the projection data.
if strcmp(tilt_type, 'tif')
    projs = ReadTIFF(projs_name);
    
elseif strcmp(tilt_type, 'mrc')
    % Permute the second and third dimensions for CS-ET convenience.
    projs = double(ReadMRC(projs_name));
    projs = permute(projs, [1, 3, 2]);
    
elseif strcmp(tilt_type, 'mat')
    % Awkward construction. Thanks MATLAB!
    file_in = load(projs_name);
    file_fields = fieldnames(file_in);
    field = file_fields{1};
    projs = getfield(file_in, field); %#ok<GFLD>
    
else
    % Something went wrong - we shouldn't be able to get here.
    error('Invalid tilt series file name. How did you trigger this error?');
end

% Load the recdata. Saved awkwardly as a field in another struct, so extract it.
var = load(recdata_name);
var_fields = fieldnames(var);
field = var_fields{1};
recdata = getfield(var, field); %#ok<GFLD>
end
