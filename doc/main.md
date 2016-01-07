# [Documentation](#top)
The CS-ET library consists of three modules:

+ The [**core**](#core) module contains the minimum set of functions needed to create 3D CS-ET tomographic reconstructions from (preprocessed) STEM tilt series, as well as 3D WBP reconstructions used for comparison.

+ The [**analysis**](#analysis) module contains functions to analyze and compare tomograms.

+ The [**im-tools**](#im-tools) module contains functions to view tomograms, as well as to create images and animations from tomograms.

This repository also contains several [**example**](#example) scripts, which demonstrate CS-ET and WBP reconstructions of the datasets used in (Guay et al., 2016).

Additionally, this repository contains several sample datasets, used by the example scripts and also useful for custom tests of the CS-ET algorithm.

---

## [](#core)[Core](core.md#top)

#### [cset](core.md#cset)
**cset** creates a 2D image array reconstructed from Radon transform data (i.e. parallel-beam STEM tilt series data) via CS-ET - an *L<sup>1</sup>*-regularized least-squares reconstruction algorithm.

#### [cset3](core.md#cset3)
**cset3** performs a CS-ET reconstruction of a 3D tomographic volume from a STEM tilt series, by decomposing the 3D volume into parallel 2D slices normal to the tilt series axis of rotation.

#### [cset\_parameters](core.md#cset_parameters)
**cset\_parameters** creates a struct containing all of the parameters needed to perform a 2D or 3D CS-ET reconstruction.

#### [get\_projs](core.md#get_projs)
**get\_projs** reads in a TIF, MRC, or MAT file containing preprocessed STEM tilt series data, along with an additional struct containing information necessary for producing a tomographic reconstruction.

#### [radon3](core.md#radon3)
**radon3** Computes the Radon transform (projections) of a 3D dataset for a given set of projection angles, using the [ASTRA Toolbox](https://github.com/astra-toolbox/astra-toolbox/). Supports CPU as well as CUDA GPU computation.

#### [setup](core.md#setup)
**setup** augments the MATLAB search path to include all directories needed by the CS-ET library.

#### [setup\_pool](core.md#setup_pool)
**setup\_pool** is a convenience wrapper for MATLAB's parpool() function, used to simplify working with parallelized and non-parallelized runs of the CS-ET algorithm.

#### [wbp3](core.md#wbp3)
**wbp3** performs a weighted backprojection (WBP) reconstruction of a 3D tomographic volume from parallel-beam STEM tilt series measurements. Uses MATLAB's Image Processing Toolbox.

---

## [](#analysis)[Analysis](analysis.md#top)

#### [align\_tomogram\_subvolume](analysis.md#align_tomogram_subvolume)
**align\_tomogram\_subvolume** uses cross-correlation to find a subvolume of one tomogram that optimally matches another, smaller tomogram.

#### [all\_compressibilities](analysis.md#all_compressibilities)
**all\_compressibilities** computes the relative compressibility of each 2D slice
(third coordinate) of an input 3D array in the identity, TV, and DB8
wavelet domains at a set of specified thresholds.

#### [compute\_compressibility](analysis.md#compute_compressibility)
**compute\_compressibility** computes the relative compressibility of an input array for a set of specified thresholds.

#### [tv](analysis.md#tv)
**tv** computes TV(u), the anisotropic total variation of an input 2D or 3D array u.

---

## [](#im-tools)[Im-tools](im-tools.md#top)

#### [add\_gaussian\_noise](im-tools.md#add_gaussian_noise)
**add\_gaussian\_noise** adds Gaussian noise to an array, with a specified signal-to-noise ratio.

#### [create\_xy\_patches](im-tools.md#create_xy_patches)
**create\_xy\_patches** creates one or more whitened x-y image patches from a
collection of tomograms. This is used for creating a visual comparison of
the outputs of multiple reconstruction methods from the same data.

#### [create\_xz\_patches](im-tools.md#create_xz_patches)
**create\_xz\_patches** creates one or more whitened x-z image patches from a collection of tomograms. This is used for creating a visual comparison of the outputs of multiple reconstruction methods from the same data.

#### [imprep](im-tools.md#imprep)
**imprep** prepares an array for use with imwrite() by rescaling the array range to [0,1]. Array values can also be clamped to a specified range before rescaling.

#### [rec\_to\_tif](im-tools.md#rec_to_tif)
**rec\_to\_tif** saves a tomogram reconstruction volume as a TIF stack.

#### [snapshot](im-tools.md#snapshot)
**snapshot** saves an image of one 2D x-y slice of a 3D input, or
the entirety of a 2D input, to a PNG file.

#### [stack\_to\_gif](im-tools.md#stack_to_gif)
**stack\_to\_gif** saves an [M,N,P] array as an [M,N] GIF with P or 2P-1 frames, depending on loop behavior.

#### [sv](im-tools.md#sv)
**sv** (**s**lice **v**iewer) allows the user to scroll back and forth through a 3D volume to view 2D slices, using horizontal motions of the mouse.

#### [whiten](im-tools.md#whiten)
**whiten** transforms an input array to have zero mean and unit variance.

---

## [](#example)[Example](example.md#top)

#### [example1\_phantom\_nano](example.md#example1_phantom_nano)
**example1\_phantom\_nano** creates CS-ET and WBP reconstructions of the simulated nanoparticle phantom used in (Guay et al., 2016).

#### [example2\_phantom\_simple](example.md#example2_phantom_simple)
**example2\_phantom\_simple** creates CS-ET and WBP reconstructions of the simple membrane phantom used in (Guay et al., 2016).

#### [example3\_phantom\_complex](example.md#example3_phantom_complex)
**example3\_phantom\_complex** creates CS-ET and WBP reconstructions of the complex membrane phantom used in (Guay et al., 2016).

#### [example4\_brightfield](example.md#example4_brightfield)
**example4\_brightfield** creates CS-ET and WBP reconstructions of the experimental bright-field dataset used in (Guay et al., 2016).

#### [example5\_darkfield](example.md#example5_darkfield)
**example5\_darkfield** creates CS-ET and WBP reconstructions of the experimental dark-field dataset used in (Guay et al., 2016).

---

[**Top**](#top)

#### [Readme](../README.md)

#### [Home](https://github.com/heyitsguay/cset)
