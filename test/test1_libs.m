function pass = test1_libs()
% TEST1_LIBS verifies that all libraries required by the CS-ET library can
% be found and are working correctly.
%
% Created: 12/31/2015
% =======
%
% Modified: 12/31/2015 "Created."
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
% pass = TEST1_LIBS() runs a series of tests to verify that all libraries
% required by the CS-ET project can be found in the search path and are
% working correctly. Returns true if all tests pass, false otherwise.
%
% Input: None
% =====
%
% Output:
% ======
% pass - If all tests pass, returns true. Otherwise, false.

tic

% Count the number of tests passed.
passes = 0; 

% Counts the number of tests failed.
fails = 0; 

% Check for a CUDA-capable GPU.
can_use_gpu = (gpuDeviceCount > 0);

% Basic setup.
try
    setup();
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('Function setup() failed with error:');
    disp(getReport(me));
end

% Test the ASTRA Toolbox, using sample scripts copied from the ASTRA
% Toolbox's sample/matlab directory.
% Save root CS-ET directory.
root_dir = cd;
try
    cd('test/libs');
    run s004_cpu_reconstruction
    % Run an additional test if CUDA-capable GPUs are detected
    if can_use_gpu
        run s003_gpu_reconstruction
    end
    drawnow;
    cd(root_dir);
    passes = passes + 1;
catch me
    fails = fails + 1;
    cd(root_dir);
    disp('ASTRA Toolbox tests failed with error:');
    disp(getReport(me));
end

% Test the Rice Wavelet Toolbox (RWT). This basically just makes sure
% MATLAB can find it. See the RWT documentation if you wish to run more
% comprehensive tests included with that toolbox.
try
    % Little test signal.
    x = makesig('LinChirp', 8);
    % Compute DB4 wavelet scaling filter
    h = daubcqf(4, 'min');
    % Number of wavelet levels.
    L = 2;  % For 8 values in x we would normally be L=
    % Compute DWT.
    [y, ~] = mdwt(x, h, L);
    % The correct values.
    y_corr = [1.1097 0.8767 0.8204 -0.5201 -0.0339 0.1001 0.2201 -0.1401];
    % Make sure residual between y and y_corr is small.
    if norm(y - y_corr) > 1e-4
        error('RWT wavelet transform is not functioning properly.');
    end
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('Rice Wavelet Toolbox tests failed with error:');
    disp(getReport(me));
end

% Test the EMIO library, used for reading MRC (and TIF) data files.
try
    ReadMRC('test/libs/test.mrc');
    ReadTIFF('test/libs/test.tif');
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('Library emio tests failed with error:');
    disp(getReport(me));
end

% Test the rot90_3D library, used for rotating 3D arrays.
try
    x = zeros(2, 3, 4);
    rot90_3D(x, 2, 3);
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('Library rot90_3D test failed with error:');
    disp(getReport(me));
end

elapsed = toc;

% Print basic GPU capabilities.
if can_use_gpu
    disp('CUDA-capable GPU detected. GPU capabilities available.');
else
    disp('No CUDA-capable GPU detected. GPU capabilities are not available.');
end

% Print test results.
fprintf(['test1_libs complete in %g s.\n' ...
         '%u tests passed, %u tests failed.\n\n'], ...
         elapsed, passes, fails);

% Close the 3 figures that the tests created.
close, close, close;

% Output whether all tests passed or not.
pass = (fails == 0);
end
