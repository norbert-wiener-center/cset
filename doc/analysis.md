# [Analysis](#top)
*A set of functions for the analysis and comparison of tomographic reconstructions.*

#### [Back](main.md#top)

## [Contents](#contents)
[align\_tomogram\_subvolume](#align_tomogram_subvolume)

[all\_compressibilities](#all_compressibilities)

[compute\_compressibility](#compute_compressibility)

[tv](#tv)

---

## [align\_tomogram\_subvolume](#align_tomogram_subvolume)
**align\_tomogram\_subvolume** uses cross-correlation to find a subvolume of one tomogram that optimally matches another, smaller tomogram.

#### Usage
subvolume = **align\_tomogram\_subvolume**(tom\_large, tom\_small) returns a subvolume of tom\_large with the same dimensions as tom\_small, which maximizes the correlation between the whitened subvolume and the whitened tom\_small. The returned (subvolume) is not whitened.

[subvolume, xcorrs] = **align\_tomogram\_subvolume**(tom\_large, tom\_small) does the same as the previous usage, but also returns a vector (xcorrs) of cross-correlations between the whitened subvolumes and the whitened tom\_small.

#### Input
Variable | Description
--- | ---
tom\_large | The larger tomogram, size [M1,N,P] where M1 > M2.             **align\_tomogram\_subvolume** returns a subvolume of tom\_large.
tom\_small | The smaller tomogram, size [M2,N,P] where M2 < M1.

#### Output
Variable | Description
--- | ---
subvolume | [M2,N,P] subvolume of tom\_large which, when whitened, maximizes the correlation with the whitened tom\_small out of all subvolumes of tom\_large with the same dimensions.
xcorrs    | [M1-M2+1,1] vector of cross-correlations between all [M2,N,P] whitened subvolumes of tom\_large and the whitened tom\_small.

---


## [all\_compressibilities](#all_compressibilities)
**all\_compressibilities** computes the relative compressibility of each 2D x-z slice
of an input 3D array in the identity, TV, and DB8 wavelet domains at a set of specified thresholds.

#### Usage
[icomp, tvcomp, wcomp] = **all\_compressibilities**(data, offset, thresholds) computes the compressibility of each 2D slice (third coordinate) of input 3D array (data) in the wavelet, TV, and identity domains at the thresholds specified in input vector (thresholds). The input (offset) is used to calculate the proper compressibility in the identity domain, in the case that it has a non-zero background value.

#### Input
Variable | Description
--- | ---
data       | 3D input array. This function calculates the 2D             compressibility of each of the P slices data(:, :, 1), ...,              data(:, :, P) in the identity, TV, and DB8 wavelet domains.
offset     | Scalar offset value, should be equal to the background value             of the input (data).
thresholds | [T,1] vector of threshold values used in the compressibility             calculations. Should be a proportion between 0 and 1.

#### Output
Variable | Description
--- | ---
icomp  | [T,P] array of compressibility values in the identity domain.
tvcomp | [T,P] array of compressibility values in the TV domain.
wcomp  | [T,P] array of compressibility values in the DB8 wavelet         domain.

---


## [compute\_compressibility](#compute_compressibility)
**compute\_compressibility** computes the relative compressibility of an input array for a set of specified thresholds.

#### Usage
com_vals = **compute\_compressibility**(u, thresholds) computes the relative
compressibility of array u for each of the relative thresholds specified
in (thresholds). For a given threshold t, this returns the proportion of elements of u with magnitude less than t*max(abs(u(:))).

#### Input
Variable | Description
--- | ---
u          | Input array to calculate the compressibility of.
thresholds | [T,1] vector of threshold values, specified as a proportion in the range (0,1).

#### Output
Variable | Description
--- | ---
com_vals | [T,1] vector of compressibility values.

---


## [tv](#tv)
**tv** computes TV(u), the anisotropic total variation of input 2D or 3D array (u). For 2D u, if *D<sub>x</sub>(u)* is the horizontal forward difference operator and  *D<sub>y</sub>(u)* is the vertical forward difference operator, then

![Definition of TV(u)](tveq.png)

The definition for 3D u is analogous, with an additional forward difference term *D<sub>z</sub>(u)* in the third dimension.

#### Usage
tvu = **tv**(u) returns TV(u), the anisotropic total variation of u.

#### Input
Variable | Description
--- | ---
u | 2D or 3D array.

#### Output
Variable | Description
--- | ---
tvu | TV(u), the anisotropic total variation of u.

---

[**Top**](#top)

#### [Back](main.md#top)
