# [Im-tools](#top)
*A collection of functions for creating and viewing images and animations of tomogram volumes.*

#### [Back](main.md#top)

## [Contents](#contents)
[add\_gaussian\_noise](#add_gaussian_noise)

[create\_xy\_patches](#create_xy_patches)

[create\_xz\_patches](#create_xz_patches)

[imprep](#imprep)

[rec\_to\_tif](#rec_to_tif)

[snapshot](#snapshot)

[stack\_to\_gif](#stack_to_gif)

[sv](#sv)

[whiten](#whiten)

---

## [add\_gaussian\_noise](#add_gaussian_noise)
**add\_gaussian\_noise** adds Gaussian noise to an array, with a specified signal-to-noise ratio.

#### Usage
xn = **add\_gaussian\_noise**(x, snr) adds Gaussian noise with variance equal to rms(x)^2 / (snr) to each component of the input array (x), where rms(x) is the root mean square of x.

xn = **add\_gaussian\_noise**(x, snr, seed) does the same as the previous usage, but seeds MATLAB's RNG with input (seed), so that results may be easily reproduced.

#### Input
Variable | Description
--- | ---
x   | Array to which the noise is added.
snr | Desired signal-to-noise ratio.
seed | (OPTIONAL) RNG seed.

#### Output
Variable | Description
--- | ---
xn | Noisy output array.

---


## [create\_xy\_patches](#create_xy_patches)
**create\_xy\_patches** creates one or more whitened x-y image patches from a
collection of tomograms. This is used for creating a visual comparison of
the outputs of multiple reconstruction methods from the same data.

#### Usage
patches = **create\_xy\_patches**(x1, y1, z1, M, N, recs) creates 2D x-y
patches from each of the tomographic reconstructions in the cell array
recs. Each patch is taken from the x-y plane with z coordinate z1, and
spans indices (y1, x1):(y1 + M, x1 + N).

patches = **create\_xy\_patches**(x1, y1, z1, z2, M, N, recs) creates 3D x-y
patches from each of the tomographic reconstructions in the cell array
recs. Each patch is a 3D stack of 2D x-y image patches with z coordinates
from z1:z2. Each 2D patch spans indices (y1, x1):(y1 + M, x1 + N). Useful
for making GIFs.

#### Input
Variable | Description
--- | ---
x1   | Top-left x coordinate of the patches.
y1   | Top-left y coordinate of the patches.
z1   | Smallest (or only) z coordinate of the patches.
z2   | (OPTIONAL) Largest z coordinate of the patches.
M    | Patch height.
N    | Patch width.
recs | [nt,1] cell array containing the 3D tomographic reconstructions from which the patches are extracted.

#### Output
Variable | Description
--- | ---
patches | [nt,1] cell array containing the image patches.

---


## [create\_xz\_patches](#create_xz_patches)
**create\_xz\_patches** creates one or more whitened x-z image patches from a collection of tomograms. This is used for creating a visual comparison of the outputs of multiple reconstruction methods from the same data.

#### Usage
patches = **create\_xz\_patches**(x1, z1, y1, M, N, recs) creates 2D x-z
patches from each of the tomographic reconstructions in the cell array
recs. Each patch is taken from the x-z plane with y coordinate y1, and
spans indices (z1, x1):(z1 + M, x1 + N).

patches = **create\_xz\_patches**(x1, z1, y1, y2, M, N, recs) creates 3D x-z
patches from each of the tomographic reconstructions in the cell array
recs. Each patch is a 3D stack of 2D x-z image patches with y coordinates
from y1:y2. Each 2D patch spans indices (z1, x1):(z1 + M, x1 + N). Useful for making GIFs.

#### Input
Variable | Description
--- | ---
x1   | Top-left x coordinate of the patches.
z1   | Top-left z coordinate of the patches.
y1   | Smallest (or only) y coordinate of the patches.
y2   | (OPTIONAL) Largest y coordinate of the patches.
M    | Patch height.
N    | Patch width.
recs | Cell containing the 3D tomographic reconstructions from which the patches are extracted.

#### Output
Variable | Description
--- | ---
patches | Cell array containing the patches, one cell for each cell in          input recs.

---


## [imprep](#imprep)
**imprep** prepares an array for use with imwrite() by rescaling the array range to [0,1]. Array values can also be clamped to a specified range before rescaling.

#### Usage
y = **imprep**(x) rescales the values in x to lie in the interval [0,1].

y = **imprep**(x, lims) first clamps the values of x to the interval
[lims(1), lims(2)], then rescales the clamped array to the interval
[0,1].

#### Input
Variable | Description
--- | ---
x    | Input array. Values get rescaled to the interval [0,1].
lims | (OPTIONAL) [2,1] vector of values to clamp x's range to before  rescaling to [0,1]. 

#### Output
Variable | Description
--- | ---
y | Rescaled array, with values in the range [0,1].

---


## [rec\_to\_tif](#rec_to_tif)
**rec\_to\_tif** saves a tomogram reconstruction volume as a TIF stack.

#### Usage
**rec\_to\_tif**(rec, file\_name) saves the tomographic reconstruction (rec) to a
TIF file with name file\_name.tif. Prior to saving, the tomogram is
rotated so that each element of the 3D TIF stack corresponds to one x-y
orthogonal view of the reconstruction volume.

**rec\_to\_tif**(rec, file\_name, rot\_flag) does the same as the previous
usage, but the boolean input rot\_flag can be used to specify whether the
reconstruction volume should be rotated prior to saving. If not rotated,
each element of the 3D TIF stack corresponds to one x-z slice of the
reconstruction volume.

**rec\_to\_tif**(rec, file\_name, rot\_flag, contrast\_scale) does the same as
the previous usage, but allows the user to clamp the range of the
reconstruction volume to an interval of values specified in input
contrast\_scale.

#### Input
Variable | Description
--- | ---
rec            | 2D or 3D tomographic reconstruction volume.
file\_name      | Name of the TIF file to be created. 
rot\_flag       | (OPTIONAL=true) Specifies whether the reconstruction volume should be rotated around its x axis prior to saving.
contrast\_scale | (OPTIONAL) [2,1] vector. Specifies a range of values to which (rec)'s range should be clamped prior to saving.

#### Output
None

---


## [snapshot](#snapshot)
**snapshot** saves an image of one 2D x-y slice of a 3D input, or
the entirety of a 2D input, to a PNG file.

#### Usage
**snapshot**(file\_name, tom, idx) with ndims(tom)==3, rotates the 3D tomogram
(tom) about its x axis to obtain an x-y orthogonal view of the volume,
and saves an image of the x-y slice with z-coordinate (idx) to a PNG file
with name (file\_name). If ndims(tom)==2, **snapshot** saves the 2D array
(tom) to a PNG file with name (file\_name).

**snapshot**(file\_name, tom, idx, min\_val, max\_val) does the same as the
previous usage, but first clamps the range of the tomogram values to
[min\_val, max\_val].

#### Input
Variable | Description
--- | ---
file\_name | Name of the PNG file to which **snapshot** saves the indicated            tomogram slice.
tom       | 2D or 3D tomographic reconstruction. Could be any other data            array of the same dimensions, I guess.
idx       | If input (tom) is a 3D array, this specifies the z-coordinate            of the x-y slice to save an image of.
min\_val   | (OPTIONAL) minimum value to clamp (tom)'s range to.
max\_val   | (OPTIONAL) maximum value to clamp (tom)'s range to.

#### Output
None

---


## [stack\_to\_gif](#stack_to_gif)
**stack\_to\_gif** saves an [M,N,P] array as an [M,N] GIF with P or 2P-1 frames, depending on loop behavior.

#### Usage
**stack\_to\_gif**(file\_name, stack) converts the [M,N,P] array (stack) into an
[M,N] gif with P frames, with file name (file\_name) and a 1/15 second
delay between each frame.

**stack\_to\_gif**(file\_name, stack, delay) does the same as the previous
usage, and allows the user to specify a custom inter-frame delay with
using input (delay).

**stack\_to\_gif**(file\_name, stack, delay, final\_delay) does the same as the
previous usage, and allows the user to specify a different inter-frame
delay after the final frame, allowing the animation to pause for a bit
longer at the final frame before looping.

**stack\_to\_gif**(file\_name, stack, delay, final\_delay, add\_reversed) does the
same as the previous usage, and allows the user to set a boolean flag
(add\_reversed) indicating if the GIF should cycle from the end frame to
the first frame in reverse.

**stack\_to\_gif**(file\_name, stack, delay, final\_delay, add\_reversed, scale)
does the same as the previous usage, and allows the user to set a custom
scaling coefficient. The imwrite function requires all images ranges to
be in [0, 255], and will clamp values to that range implicitly. The
default function procedure is to scale the max value in (stack) to 255.

#### Input
Variable | Description
--- | ---
file\_name    | Name of the GIF to be created.
stack        | [M,N,P] array to be converted to a GIF.
delay        | (OPTIONAL=1/15) Delay time between GIF frames, in seconds.
final\_delay  | (OPTIONAL=3) Delay time after the final GIF frame before               looping, in seconds.
add\_reversed | (OPTIONAL=false) If true, the GIF will animate from the               last frame in (stack) to the first frame in reversed               order.
scale        | (OPTIONAL) Multiplicative scaling parameter applied to               (stack).

#### Output
None

---


## [sv](#sv)
**sv** (**s**lice **v**iewer) allows the user to scroll back and forth through a 3D
volume to view 2D slices, using horizontal motions of the mouse.

#### Usage
**sv**(vol) creates a new figure with callback functions associated with
mouse motion and key presses. Horizontal mouse motion will change which
2D slice of the 3D volume (vol) is displayed (indexed by its third
coordinate), while a key press will exit the function.

**sv**(vol, mag\_level) does the same as the previous usage, but allows the
user to specify a magnification level through input (mag\_level).
Magnification level should be specified as a percentage, e.g. the default
100% magnification corresponds to mag\_level = 100.

**sv**(vol, mag\_level, im\_range) does the same as the previous usage, but
allows the user to specify a display range when displaying each slice
using imshow(). The default value is [min(vol(:)), max(vol(:))].

#### Input
Variable | Description
--- | ---
vol       | A 3D array. **sv** displays slices of this array, indexed by its           third coordinate.
mag\_level | (OPTIONAL=100) The magnification level of the displayed            slices. Should be input as a percentage.
im\_range  | (OPTIONAL) The image display range, passed to imshow() and            operates the same as the [low high] input to the imshow()    function.

#### Output
None

---


## [whiten](#whiten)
**whiten** transforms an input array to have zero mean and unit variance.

#### Usage
w = **whiten**(x) transforms (x) to have zero mean and unit variance and outputs the transformed array as (w).

#### Input
Variable | Description
--- | ---
x | An array of arbitrary dimensions.

#### Output
Variable | Description
--- | ---
w | A transformed version of input (x) with zero mean and unit variance.

---

[**Top**](#top)

#### [Back](main.md#top)
