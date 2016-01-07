function pass = test3_analysis()
% TEST3_ANALYSIS tests the functions in the 'analysis' module of the CS-ET
% library - a set of functions used to analyze and compare tomogram
% volumes.
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
% pass = TEST3_ANALYSIS() runs a series of tests to verify that the
% functions in the 'analysis' module of this library are working correctly.
% Returns true if all tests pass, false otherwise.
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

% Test tomogram alignment function.
try
    % Load a test volume.
    x = load('test/cset/test_vol.mat');
    vol = x.vol;
    
    % Make a smaller version of vol.
    vol_small = vol(14:(end-6), :, :);
    
    % Add a bit of noise to vol_small.
    rng(1126); % Specify a RNG seed, for deterministic results.
    vol_noisy = vol_small + 0.05 * randn(size(vol_small));
    
    % Find the subvolume of vol that optimally aligns with vol_small.
    [subvolume, ~] = align_tomogram_subvolume(vol, vol_noisy);
    
    % Check that we got the right subvolume.
    err = norm(subvolume(:) - vol_small(:));
    if err > 0
        error('Incorrect subvolume selected.');
    end
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('Subvolume alignment test failed with error:');
    disp(getReport(me));
end

% Test compressibility analysis.
try
    [ic, tvc, wc] = all_compressibilities(vol, 0.1, [0, 0.01, 0.05, 0.1]);
    
    % Concatenate outputs, to simplify testing.
    c_vals = cat(3, ic, tvc, wc);
    
    % Load stored compressibility values.
    x = load('test/analysis/test_c_vals.mat');
    c_vals_correct = x.c_vals;
    
    % Make sure the compressibility calculation worked.
    assert(norm(c_vals(:)-c_vals_correct(:)) < 1e-6, ...
        'Calculated compressibility values do not match reference version.');
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('Compressibility analysis test failed with error:');
    disp(getReport(me));
end

elapsed = toc;

% Print results.
fprintf(['test3_analysis complete in %g s.\n' ...
         '%u tests passed, %u tests failed.\n\n'], ...
         elapsed, passes, fails);
     
% Output whether all tests passed or not.
pass = (fails == 0);
end