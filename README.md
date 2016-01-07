# cset
This is a MATLAB project for compressed sensing (CS) reconstruction of electron tomogram (ET) volumes. The tools in this project are useful for creating 3D tomogram volumes from parallel beam tilt series created through scanning transmission electron microscopy (STEM), though it should be readily adaptable to other forms of parallel beam tomography as well. The CS-ET algorithm is suitable for the reconstruction of sparse datasets from undersampled tilt series, as well as for creating denoised tomograms from fully-sampled tilt series. See (Guay et al., 2016), available at [*insert link here when it's available*], for more details.

## Example
```Matlab
% Adds the necessary directories to the search path, including a 
% path to the dataset in the directory '../cset-data/phantom-complex'.
setup('phantom-complex');

% Load the (simulated ET) projection data for this dataset, along
% with a struct (recdata) containing additional info needed for the 
% CS-ET reconstruction - projection angles, background intensity, 
% and reconstruction dimensions.
[projs, recdata] = get_projs();

% Split-Bregman optimization hyperparameters.
mu = 0.005;
lambda = mu * 10;
gamma = mu * 20;
kappa = mu * 0.1;
% Inner and outer (Bregman) loop counts.
n_inner = 20;
n_outer = 30;

% Conjugate gradient parameters.
cg_tol = 1e-4;
max_cgiter = 12;

% Flag to use CUDA for Radon transforms.
gpu_flag = true;

% If true, run multiple 2D reconstructions in parallel
p_flag = false;
% Number of parallel pool workers.
pool_size = 26;
% If p_flag is true, this will call MATLAB's parpool function.
setup_pool(p_flag, pool_size);

% Create a CS-ET parameter struct.
params = cset_parameters(recdata, mu, lambda, gamma, kappa, n_inner, n_outer, p_flag, cg_tol, max_cgiter, gpu_flag);

% Create a 2D CS-ET reconstruction of a single x-z slice of 
% the dataset (y=7).
rec2d = cset(projs(:, :, 7), params);

% Create a 3D CS-ET reconstruction.
rec3d = cset3(projs, params);
```

Additional, fleshed-out examples are included in the [example](example/) directory. More information about these examples can be found [here](doc/example.md).

## Installation
This project requires the MATLAB Image Processing Toolbox, along with two other third-party MATLAB libraries, the [ASTRA Toolbox](https://github.com/astra-toolbox/astra-toolbox/) and the [Rice Wavelet Toolbox](https://github.com/ricedsp/rwt). Please install those libraries as directed.  Additional small MATLAB libraries are included in this repo, in the [lib](lib/) folder.

To use optional GPU capabilities with the ASTRA Toolbox, a CUDA-capable GPU is required, along with the associated CUDA drivers. See the [CUDA Downloads](https://developer.nvidia.com/cuda-downloads) page for information on installing CUDA drivers, **prior** to the installation of the ASTRA Toolbox. Using optional multithreading capabilities requires the MATLAB Parallel Processing Toolbox.

Once that is done, this library needs to know the locations of the $ASTRA\_ROOT/matlab folder and the $RWT\_ROOT/bin folder, as well as a default location for any datasets you wish to use with this project. These are specified in the [settings.cfg](settings.cfg) file, located in the project root directory. 

## Tests
Test functions can be found in the [test](test/) directory. 

`run_tests()` runs all tests.

`test1_libs()` verifies that all necessary libraries can be found and are functioning properly. 

`test2_core()` verifies that the core CS-ET algorithms are functioning properly.

`test3_analysis()` verifies that the 'analyis' module is functioning properly.

`test4_imtools()` verifies that the 'im-tools' module is functioning properly.

## Documentation
The documentation for this library can be viewed [here](doc/main.md).

## Contributors
Contributions are welcome! Please use GitHubs Issues tracker to report any bugs you may find. 

## References
To come.

## License
Copyright (c) 2015-2016 Matthew Guay

This software is licensed under the GPLv3 License.

## Contact
Contact Matt Guay at mguay@math.umd.edu with any comments or questions you may have.
