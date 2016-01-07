function pass = test4_imtools()
% TEST4_IMTOOLS tests the functions in the 'im-tools' module of the CS-ET
% library - a set of functions to create images and animations of 3D
% arrays, i.e. tomogram volumes.
%
% Created: 01/01/2016
% =======
%
% Modified: 01/01/2016 "Created."
% ========
%
% Author: Matthew Guay
% ======  mguay@math.umd.edu
%         Applied Mathematics & Statistics, and Scientific Computation
%         Department of Mathematics
%         University of Maryland, College Park
%         Copyright (C) 2016
%
% Usage:
% =====
% pass = TEST4_IMTOOLS() runs a series of tests to verify that the
% functions in the im-tools module are working correctly. Returns true if
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

% Test create_xy_patches() and create_xz_patches().
try
    % Load test volumes
    x = load('test/cset/test_rec_cset.mat');
    rec_cset = x.rec_cset;
    x = load('test/cset/test_rec_wbp.mat');
    rec_wbp = x.rec_wbp;
    
    % Create x-y patches.
    pxy_cell = create_xy_patches(10, 1, 20, 16, 16, {rec_cset, rec_wbp});
    patches_xy = [pxy_cell{1} pxy_cell{2}];
    
    % Create x-z patches.
    pxz_cell = create_xz_patches(10, 10, 8, 16, 16, {rec_cset, rec_wbp});
    patches_xz = [pxz_cell{1} pxz_cell{2}];
    
    % Load reference patches for comparison.
    x = load('test/im-tools/test_patches.mat');
    patches_xy_correct = x.patches_xy;
    patches_xz_correct = x.patches_xz;
    
    % Make sure the patches were created correctly.
    assert(norm(patches_xy(:)-patches_xy_correct(:)) == 0, ...
        'Created x-y patches do not match reference version.');
    assert(norm(patches_xz(:)-patches_xz_correct(:)) == 0, ...
        'Created x-z patches do not match reference version.');
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('Patch creation test failed with error:');
    disp(getReport(me));
end

% Test imprep().
try
    % Take a 2D slice from rec_cset.
    img = rec_cset(:, :, 8);
    
    % Run img through imprep
    primg = imprep(img, [0.1, 0.5]);
    
    % Load reference version.
    x = load('test/im-tools/test_primg.mat');
    primg_correct = x.primg;
    
    % Compare with reference version.
    assert(norm(primg(:)-primg_correct(:)) < 1e-6, ...
        'Image created with imprep() does not match reference version.');
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('imprep() test failed with error:');
    disp(getReport(me));
end

% Test rec_to_tif().
try
    % Create a TIF file of rec_cset.
    rec_to_tif(rec_cset, 'cset.tif', false, [0.1, 0.7]);
    
    % Load cset.tif, and a reference version.
    cset_tif = ReadTIFF('cset.tif');
    cset_tif_correct = ReadTIFF('test/im-tools/cset_correct.tif');
    
    % Verify that the files match. I suppose additional hidden errors could
    % creep in if something were wrong with the ReadTIFF() function, as
    % well...
    assert(norm(cset_tif(:)-cset_tif_correct(:)) == 0, ...
        'Output of rec_to_tif() does not match reference version.');
    
    % Delete the created TIF file.
    delete('cset.tif');
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    delete('cset.tif');
    disp('rec_to_tif() test failed with error:');
    disp(getReport(me));
end

% Test snapshot().
try
    % Create a snapshot of an x-y slice of rec_cset.
    snapshot('cset', rec_cset, 8, [0.2, 0.5]);
    
    % Load cset.png and a reference version.
    cset_png = double(imread('cset.png'));
    cset_png_correct = double(imread('test/im-tools/cset_correct.png'));
    
    % Verify that the files match.
    assert(norm(cset_png(:)-cset_png_correct(:)) == 0, ...
        'Output of snapshot() does not match reference version.');
    
    % Delete the created PNG file.
    delete('cset.png');
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    delete('cset.png');
    disp('snapshot() test failed with error:');
    disp(getReport(me));
end

% Test stack_to_gif().
try
    % Create a gif from rec_cset.
    stack_to_gif('cset.gif', rec_cset, 1/15, 2, true);
    
    % Load cset.gif and a reference version.
    cset_gif = double(imread('cset.gif'));
    cset_gif_correct = double(imread('test/im-tools/cset_correct.gif'));
    
    % Verify that the files match.
    assert(norm(cset_gif(:)-cset_gif_correct(:)) == 0, ...
        'Output of stack_to_gif() does not match reference version.');
    
    % Delete the created GIF file.
    delete('cset.gif');
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    delete('cset.gif');
    disp('stack_to_gif() test failed with error:');
    disp(getReport(me));
end

% You might think we'd next test sv() here. However, due to its dependence
% on user mouse and keyboard input, and my unwillingness to figure out how
% to use the java.awt.Robot library concurrently to simulate that input,
% there will be no sv() test for now. Don't mess with it, and it should be
% fine.

% Test whiten()
try
    % Pull a slice from rec_cset to whiten.
    slice = rec_cset(:, :, 9);
    
    % Whiten the slice.
    wslice = whiten(slice);
    
    % Load the reference version.
    x = load('test/im-tools/test_wslice.mat');
    wslice_correct = x.wslice;
    
    % Verify that the slices match.
    assert(norm(wslice(:)-wslice_correct(:)) < 1e-6, ...
        'Output of whiten() does not match reference version.');
    
    passes = passes + 1;
catch me
    fails = fails + 1;
    disp('whiten() test failed with error:');
    disp(getReport(me))
end

elapsed = toc;

% Print test results.
fprintf(['test4_imtools complete in %g s.\n' ...
         '%u tests passed, %u tests failed.\n\n'], ...
         elapsed, passes, fails);
     
% Output whether all tests passed or not.
pass = (fails == 0);
end    