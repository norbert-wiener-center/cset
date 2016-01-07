function run_tests()
% RUN_TESTS runs all the tests created to verify the proper functioning of
% the CS-ET library. 
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
% RUN_TESTS() runs all the tests created to verify the proper functioning
% of the CS-ET library. Should any tests fail, RUN_TESTS initiates a
% platform-independent self-destruct sequence, engulfing your computer and
% its immediate surroundings in a hellish vortex of fire and molten
% plastic. Please stand clear.
%
% Input: None
% =====
%
% Output: None
% ======

test1_libs();

test2_core();

test3_analysis();

test4_imtools();
end