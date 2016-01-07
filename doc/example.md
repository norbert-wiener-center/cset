# [Example](#top)
*A collection of scripts exhibiting CS-ET and WBP reconstructions of the datsets used in* "Compressed Sensing Electron Tomography for Determining Biological Structure" *(Guay et al., 2016).*

#### [Back](main.md#top)

## [Contents](#contents)
[example1\_phantom\_nano](#example1_phantom_nano)

[example2\_phantom\_simple](#example2_phantom_simple)

[example3\_phantom\_complex](#example3_phantom_complex)

[example4\_brightfield](#example4_brightfield)

[example5\_darkfield](#example5_darkfield)

---


## [example1\_phantom\_nano](#example1_phantom_nano)
**example1\_phantom\_nano** creates CS-ET and WBP reconstructions of the simulated nanoparticle phantom used in (Guay et al., 2016), from projections corrupted by low levels of noise.

This phantom is 2D, and is therefore the fastest to reconstruct among the datasets used in the example scripts. This script creates an image of the 1x- and 3x-undersampled CS-ET and WBP reconstructions, in addition to computing reconstructions from 2x-undersampled projection data. The image generated may be compared with [example1\_phantom\_nano.png](../example/example1_phantom-nano/example1_phantom_nano.png).

---


## [example2\_phantom\_simple](#example2_phantom_simple)
**example2\_phantom\_simple** creates CS-ET and WBP reconstructions of the simple membrane phantom used in (Guay et al., 2016). Reconstructions are made from both noiseless and noisy projection data, to demonstrate the robustness to noise of the CS-ET algorithm.

One cell of the script creates 2D CS-ET reconstructions of a single x-z slice of the simple membrane phantom, and creates an image of the reconstructions from noisy and noiseless projections. This image may be compared with [example2\_phantom\_simple.png](../example/example2_phantom-simple/example2_phantom_simple.png).

Additional cells create full 3D CS-ET and WBP reconstructions of the simple membrane phantom at various levels of undersampling. **Note**: these cells may take a long time to run without GPU and multithreading support.

---

## [example3\_phantom\_complex](#example3_phantom_complex)
**example3\_phantom\_complex** creates CS-ET and WBP reconstructions of the complex membrane phantom used in (Guay et al., 2016). Reconstructions are made from both noiseless and noisy projection data, to demonstrate the robustness to noise of the CS-ET algorithm.

One cell of the script creates 2D CS-ET reconstructions of a single x-z slice of the complex membrane phantom, and creates an image of the reconstructions from noiseless and noisy projection data. This image may be compared with [example3\_phantom\_complex.png](../example/example3_phantom-complex/example3_phantom_complex.png).

Additional cells create full 3D CS-ET and WBP reconstructions of the complex membrane phantom at various levels of undersampling. **Note**: these cells may take a long time to run without GPU and multithreading support.

---

## [example4\_brightfield](#example4_brightfield)
**example4\_brightfield** creates CS-ET and WBP reconstructions of the experimental bright-field dataset used in (Guay et al., 2016). This dataset consists of a single-axis STEM tilt series taken in the bright-field imaging mode, preprocessed using [IMOD](http://bio3d.colorado.edu/imod/).

**NOTE**: This example script requires a dataset not included in this repository. You can download it from [**INSERT LINK HERE**]

One cell of the script creates a 2D CS-ET reconstruction of a single x-z slice of the tomogram, and creates an image of this reconstruction. This image may be compared with [example4\_brightfield.png](../example/example4_brightfield/example4_brightfield.png).

An additional cell creates 3D CS-ET and WBP reconstructions of the full dataset, at various levels of undersampling. **Note**: This will take a really long time to run without ample GPU and multithreading support. As a reference, it takes about an hour to run using 26 parallel pool workers on a workstation with dual eight-core Intel(R) Xeon(R) processors and an NVIDIA Tesla(TM) K20C GPU.

---

## [example5\_darkfield](#example5_darkfield)
**example5\_darkfield** creates CS-ET and WBP reconstructions of the experimental dark-field dataset used in (Guay et al., 2016). This dataset consists of a single-axis STEM tilt series taken in the dark-field imaging mode, preprocessed using [IMOD](http://bio3d.colorado.edu/imod/).

**NOTE**: This example script requires a dataset not included in this repository. You can download it from [**INSERT LINK HERE**]

One cell of the script creates a 2D CS-ET reconstruction of a single x-z slice of the tomogram, and creates an image of this reconstruction. This image may be compared with [example5\_darkfield.png](../example/example5_darkfield/example5_darkfield.png).

An additional cell creates 3D CS-ET and WBP reconstructions of the full dataset, at various levels of undersampling. **Note**: This will take a really long time to run without ample GPU and multithreading support. As a reference, it takes about an hour to run using 26 parallel pool workers on a workstation with dual eight-core Intel(R) Xeon(R) processors and an NVIDIA Tesla(TM) K20C GPU.

---

[**Top**](#top)

#### [Back](main.md#top)
