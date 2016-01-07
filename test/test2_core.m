function pass = test2_core()
% TEST2_CORE tests the core functions of the CS-ET library - the minimum
% set needed to create 3D CS-ET (and WBP) reconstructions.
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
% pass = TEST2_CORE() runs a series of tests to verify that the core
% functionality of the CS-ET project is working correctly. Returns true if
% all tests pass, false otherwise.
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

% Basic setup.
try
    setup();
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('Function setup() failed with error:');
    disp(getReport(me));
end

% Compute Radon transforms of a little dataset.
try
    % Load a test volume.
    x = load('test/cset/test_vol.mat');
    vol = x.vol;
    
    % Load the stored Radon transform of that volume.
    x = load('test/cset/test_rad_vol.mat');
    rad_vol_correct = x.rad_vol;
    
    % Projection angles.
    theta = -50:5:50;
    
    % Perform a Radon transform.
    rad_vol = radon3(vol, theta);
    
    % Check that the Radon transform did what it's supposed to.
    assert(norm(rad_vol(:)-rad_vol_correct(:)) < 1e-6, ...
        'Calculated Radon transform does not match reference version.');
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('Radon transform calculation test failed with error:');
    disp(getReport(me));
end

% Load a tilt series and perform CS-ET and WBP reconstructions.
try
    % Load in a test tilt series and reconstruction data.
    [projs, recdata] = get_projs('test/cset/test_tilt.mat', 'test/cset/test_recdata.mat');
    
    % Split-Bregman parameters.
    % Regularization hyperparameters.
    mu = 0.005;
    lambda = mu * 10;
    gamma = mu * 20;
    kappa = mu * 1;
    % Inner and outer loop counts. Not enough for good convergence, but
    % enough to make sure the algorithm works.
    n_inner = 4;
    n_outer = 3;
    
    % Create CS-ET parameter sets.
    params = cset_parameters(recdata, mu, lambda, gamma, kappa, n_inner, n_outer);
    
    % Run a CS-ET reconstruction.
    rec_cset = cset3(projs, params);
    
    % Run a WBP reconstruction.
    rec_wbp = wbp3(projs, params);
    
    % Load stored comparision CS-ET and WBP reconstructions.
    x = load('test/cset/test_rec_cset.mat');
    rec_cset_correct = x.rec_cset;
    x = load('test/cset/test_rec_wbp.mat');
    rec_wbp_correct = x.rec_wbp;
    
    % Check that the reconstructions did what they're supposed to.
    assert(norm(rec_cset(:)-rec_cset_correct(:)) < 1, ...
        'Calculated CS-ET reconstruction does not match reference version.');
    assert(norm(rec_wbp(:)-rec_wbp_correct(:)) < 1, ...
        'Calculated WBP reconstruction does not match reference version.');
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('Tomogram reconstruction tests failed with error:');
    disp(getReport(me));
end

% Perform a parallelized CS-ET reconstruction.
try
    % First, make this false to test setup_pool() in one configuration.
    parallel_flag = false;
    % Number of parallel pool workers.
    pool_size = 2;
    % If parallel_flag is true, this will call MATLAB's parpool function.
    setup_pool(parallel_flag, pool_size);
    
    % Now make parallel_flag true to test setup_pool() in the other
    % configuration.
    parallel_flag = true;
    setup_pool(parallel_flag, pool_size);

    % New parameter set with parallel_flag.
    params2 = cset_parameters(recdata, mu, lambda, gamma, kappa, n_inner, n_outer, parallel_flag);
    
    % Run a parallelized CS-ET reconstruction.
    rec_cset = cset3(projs, params2);
    
    % Delete the parallel pool
    delete(gcp('nocreate'));
    
    % Check that the reconstruction worked correctly.
    assert(norm(rec_cset(:)-rec_cset_correct(:)) < 1, ...
        'Calculated parallel CS-ET reconstruction does not match reference version.');
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    delete(gcp('nocreate'));
    disp('Parallelized CS-ET reconstruction test failed with error:');
    disp(getReport(me));
    disp('Inessential, but you will be unable to run parallelized CS-ET reconstructions.');
end

elapsed = toc;

% Print test results.
fprintf(['test2_core complete in %g s.\n' ...
         '%u tests passed, %u tests failed.\n\n'], ...
         elapsed, passes, fails);
     
% Output whether all tests passed or not.
pass = (fails == 0);
end
