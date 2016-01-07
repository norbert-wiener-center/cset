# [Core](#top)
*The core functions of the CS-ET library.*

#### [Back](main.md#top)

## [Contents](#contents)
[cset](#cset)

[cset3](#cset3)

[cset\_parameters](#cset_parameters)

[get\_projs](#get_projs)

[radon3](#radon3)

[setup\_pool](#setup_pool)

[wbp3](#wbp3)

---


## [cset](#cset)
**cset** creates a 2D image array reconstructed from Radon transform data (i.e. parallel-beam STEM tilt series data) via CS-ET - an *L<sup>1</sup>*-regularized least-squares reconstruction algorithm. STEM data must be preprocessed using e.g. [IMOD](http://bio3d.colorado.edu/imod/) prior to use with this function.

#### Usage
u = **cset**(proj, params) performs a reconstruction using several parameters saved in a struct (params), described in the **Input** section below.

#### Input
Variable | Description
--- | ---
proj | [N,T] array of T 2D Radon transform measurements, each size N.
params  | Struct containing fields .b, .theta, .M, .N, .mu, .lambda, .gamma, .kappa, .n\_inner, .n\_outer, and optionally, .cg\_tol, .max\_cgiter, .gpu\_flag, .im\_flag, .anim\_flag, whose uses are described below.
.b | Reconstruction background value. Estimate this a priori.
.theta | [T,1] vector of Radon projection angles in the range [-90<sup>o</sup>, 90<sup>o</sup>].
.M | Reconstruction array depth (z dimension) (number of rows).
.N | Reconstruction array width (x dimension) (number of columns).
.mu | Split-Bregman parameter - controls u's update step size.
.lambda | Split-Bregman parameter - TV regularization weight.
.gamma | Split-Bregman parameter - Identity L1 regularization weight.
.kappa | Split-Bregman parameter - DWT L1 regularization weight.
.n\_inner | Number of inner loop iterations.
.n\_outer | Number of outer (Bregman) loop iterations.
.cg\_tol  | (OPTIONAL=1e-4) Conjugate gradient stopping tolerance.
.max\_cgiter | (OPTIONAL=12) Maximum CG iterations before termination.
.gpu\_flag | (OPTIONAL=true) If true, use CUDA Radon and adjoint Radon transforms. **cset** will only use GPU capabilities if MATLAB detects a CUDA-capable device.
.im\_flag | (OPTIONAL=false) If true, display u after each update iteration.
.anim\_flag | (OPTIONAL=false) If true, return an [M,N,n\_inner*n\_outer] array containing each iteration of u.

#### Output
Variable | Description 
--- | ---
u | If anim\_flag is false, an [M,N] array containing the reconstruction data at the final iteration of the algorithm. If anim\_flag is true, an [M,N,n\_inner*n\_outer] array containing the reconstruction data after every iteration of the algorithm.

---


## [cset3](#cset3)
**cset3** performs a CS-ET reconstruction of a 3D tomographic volume from a STEM tilt series, by decomposing the 3D volume into parallel 2D slices normal to the tilt series axis of rotation and running [cset](#cset) on each of them.

#### Usage
**cset3**(projs, params) performs a reconstruction from measurements (projs) using parameters saved in struct 'params'. See the **Input** section of [cset](#cset) for a definition of the fields which params must contain. An additional optional field in params is 'parallel\_flag', defined in the **Input** section below.

#### Input
Variable | Description
--- | ---
projs | [N,T,P] array of Radon transform data. P slices, each slice containing T projections of length N.
params | Struct containing fields defined in the **Input** section of [cset](#cset), and the optional field parallel\_flag described below.
.parallel_flag | Boolean value indicating if the CS-ET algorithm should perform the 2D slice reconstructions in parallel using MATLAB's parfor structure.

#### Output
Variable | Description
--- | ---
recs | [M,N,P] tomographic reconstruction volume. In agreement with the coordinate system typically used in electron microscopy, M is the reconstruction depth (z dimension), N is the reconstruction width (x dimension), P is the reconstruction length (y dimension).

---


## [cset\_parameters](#cset_parameters)
**cset\_parameters** creates a struct containing all of the parameters needed to perform a 2D or 3D CS-ET reconstruction.

#### Usage
params = **cset\_parameters**(b, theta, M, N, mu, lambda, gamma, kappa,
n\_inner, n\_outer, parallel\_flag, cg\_tol, max\_cgiter, gpu\_flag, im\_flag, anim\_flag) creates struct (params) with the first ten required fields ('b', ..., 'n\_outer'), as well as any of the optional fields. Field names are the same as the names of this function's inputs.
'
params = **cset\_parameters**(recdata, mu, lambda, gamma, kappa, n\_inner, n\_outer, parallel\_flag, cg\_tol, max\_cgiter, gpu\_flag, im\_flag, anim\_flag) creates the same output struct, but uses the input struct (recdata), which already has fields .b, .theta, .M, and .N, to create those same fields in (params). All other inputs are identical to the previous usage.

#### Input
Variable | Description
--- | ---
b             | Scalar reconstruction background value.
theta         | [T,1] vector of measurement angles.
M             | Reconstruction depth (z dimension, 1st index).
N             | Reconstruction width (x dimension, 2nd index).
mu            | Split-Bregman step-size parameter.
lambda        | Split-Bregman TV regularization hyperparameter.
gamma         | Split-Bregman Identity L1 regularization hyperparameter.
kappa         | Split-Bregman L1(DWT) regularization hyperparameter
n\_inner       | Split-Bregman number of inner loop iterations.
n\_outer       | Split-Bregman number of outer loop iterations.
parallel\_flag | (OPTIONAL=false) If true, run multiple 2D reconstructions in parallel using parfor.
cg\_tol        | (OPTIONAL=1e-4) Conjugate gradient stopping tolerance.
max\_cgiter    | (OPTIONAL=12) Conjugate gradient max iterations before termination.
gpu\_flag      | (OPTIONAL=true) If true, use CUDA Radon and adjoint Radon transforms.
im\_flag       | (OPTIONAL=false) If true, display CS-ET reconstruction after each update iteration.
anim\_flag     | (OPTIONAL=false) If true, [cset](#cset) returns an array containing each iteration of the reconstruction. Useful for creating animations.

#### Output
Variable | Description
--- | ---
params           | Struct containing a field for each of the inputs, each with the same name as the corresponding input and containing the same data.

---


## [get\_projs](#get_projs)
**get\_projs** reads in a TIF, MRC, or MAT file containing preprocessed STEM tilt series data, along with an additional struct containing information necessary for producing a tomographic reconstruction.

#### Usage
[projs, recdata] = **get\_projs**(projs\_name, recdata\_name) reads in TIF, MRC or MAT file named (projs\_name), along with additional data from MAT file (recdata\_name). Output (projs) is an [N,T,P] array projection data (P projections, each with T measurements each of length N). Output (recdata) is a struct containing five fields: .theta, .M, .N, .P, .b.

[projs, recdata] = **get\_projs**() reads in TIF, MRC or MAT files with the default names: 'tilt.mrc'/'tilt.tif'/'tilt.mat' and 'recdata.mat'. This should work if you've already called [setup](#setup) with arguments to add a data directory to the search path.

#### Input
Variable | Description
--- | ---
projs\_name   | (OPTIONAL, default='tilt.***') Name of the TIF, MRC, or MAT tilt series file to read in. If none is provided, the default is automatically determines the appropriate file extension.
recdata\_name | (OPTIONAL, default='recdata.mat') Name of the MAT file containing the recdata struct, which holds other data required for performing a reconstruction. See (recdata)'s description below in **Output** for an explanation of the fields.

#### Output
Variable | Description
--- | ---
projs   | [N,T,P] array of projection data (P projections, each with T   measurements each of length N).
recdata | Struct containing the five fields .theta, .M, .N, .P, .b, described below:
 .theta | Tx1 vector of projection angles. 
 .M | Reconstruction depth (z dimension). 
 .N | Reconstruction width (x dimension). 
 .P | Reconstruction length (y dimension). 
 .b | Reconstruction background value.
 
---
 
              
## [radon3](#radon3)
**radon3** Computes the Radon transform (projections) of a 3D dataset for a given set of projection angles, using the [ASTRA Toolbox](https://github.com/astra-toolbox/astra-toolbox/). Supports CPU as well as CUDA GPU computation.

#### Usage
projs = **radon3**(data, theta) computes the Radon transform of each 2D slice of input (data) along the third dimension of (data), at each of the angles specified in theta.

projs = **radon3**(data, theta, num\_detectors) does the same, but specifies the number of projection detectors (i.e. the width of the Radon projections) in num\_detectors.

projs = **radon3**(data, theta, num\_detectors, gpu\_flag) does the same, but specifies whether to use ASTRA's CUDA tools to perform the computation by passing the boolean input gpu\_flag.

projs = **radon3**(data, theta, num\_detectors, gpu\_flag, ss\_level) does the same, but specifies the level of supersampling to use for ASTRA's CUDA-driven Radon transform.

#### Input
Variable | Description
--- | ---
data          | [M,N,P] array containing a volume to be projected into a tilt series.
theta         | [T,1] vector of Radon transform projection 
num\_detectors | (OPTIONAL=size(data,2)) Number of detectors for each projection (i.e. width of each projection, in pixels). 
gpu\_flag      | (OPTIONAL=true) If true, use the ASTRA CUDA tools to compute Radon transforms, if MATLAB detects a CUDA-capable device.
ss\_level      | (OPTIONAL=2) Level of supersampling for the Radon transform. Only does something if gpu_flag is true.

#### Output
Variable | Description
--- | ---
tilt\_series | [num\_detectors,T,P] array. [num\_detectors,T] Radon transform data for each of the P dimension-3 slices of (data).
      
---      
      
              
## [setup](#setup)
**setup** augments the MATLAB search path to include all directories needed by the CS-ET library.

#### Usage
**setup**() performs setup without adding any data folders to the search
path.

**setup**(data\_name) performs setup and adds folder (data\_name)
located in the default data directory, which can be set in the *settings.cfg* file.

**setup**(data\_name, data\_location) performs setup and adds folder (data\_name) located in directory (data\_location) to the search path.

#### Input
Variable | Description
--- | ---
data\_name     | (OPTIONAL) String specifying the name of the dataset to use. This name should be a folder in your default data directory, or in the directory specified with input (data\_location).
data\_location | (OPTIONAL) String specifying the directory in which the folder (data\_name) is located. Default value is specified in the file *settings.cfg*.

#### Output
None

----


## [setup\_pool](#setup_pool)
**setup\_pool** is a convenience wrapper for MATLAB's parpool() function, used to simplify working with parallelized and non-parallelized runs of the [cset3](#cset3) and [wbp3](#wbp3).

#### Usage
**setup\_pool**(parallel\_flag, n) creates a MATLAB parallel pool with (n) workers if (parallel\_flag) is true and no current pool with (n) workers exists.

#### Input
Variable | Description
--- | ---
parallel\_flag | Boolean variable. If true, **setup\_pool** will create a parallel pool if none exists with (n) workers.
n             | Number of workers in the created parallel pool.

#### Output
None

----


## [wbp3](#wbp3)
**wbp3** performs a weighted backprojection (WBP) reconstruction of a 3D tomographic volume from parallel-beam STEM tilt series measurements. Uses MATLAB's Image Processing Toolbox.

#### Usage
**wbp3**(projs, params) performs a reconstruction from measurements (projs) using parameters saved in struct (params). See the **Input** section of [cset](#cset) for a definition of the fields which run must contain. An additional optional field in run is (parallel\_flag), defined in the [cset3](#cset3) **Input** section.

#### Input
Variable | Description
--- | ---
projs | [N,T,P] array of Radon transform data. P slices, each slice containing T projections of length N.
params   | Struct containing information to set up the tomogram volume reconstruction.

#### Output
Variable | Description
--- | ---
recs | [M,N,P] tomographic reconstruction volume. M is the reconstruction depth (z dimension), N is the reconstruction width (x dimension), P is the reconstruction length (y dimension).

---

[**Top**](#top)

#### [Back](main.md#top)
