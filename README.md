# cu-acetate-spin-density-fit
MATLAB code to fit a real-space spin-density model (anisotropic Gaussians) to 3D inelastic neutron scattering data of the triplet excitation in copper acetate. Accompanies arXiv:2604.13751

## What it does

Copper acetate consists of isolated Cu(II) dimers with a singlet ground state and a triplet excitation. The script:

1. loads the background-subtracted intensity of the triplet excitation, measured on a single crystal by neutron time-of-flight and stored on a 3D grid in reciprocal space;
2. excludes regions that cannot be used (empty voxels, direct beam, spurions, outer region);
3. fits the data with the Fourier transform of a spin-density model (weighted least squares, `lsqnonlin`);
4. prints the fitted parameters with uncertainties and reduced chi-squared;
5. plots data, best-fit model and residuals in the K = 0 plane and along three 1D cuts.

## The model

- The unit cell contains two differently oriented dimers (Cu1-Cu2 and Cu3-Cu4). Each Cu carries a cluster of five atoms: the Cu and its four nearest oxygens. Atomic positions come from single-crystal x-ray diffraction.
- Each atom carries a normalised, axially symmetric 3D Gaussian spin density. Its Fourier transform is

  f(Q) = exp( -1/2 ( &sigma;<sub>&perp;</sub>&sup2; Q<sub>&perp;</sub>&sup2; + &sigma;<sub>&#8741;</sub>&sup2; Q<sub>&#8741;</sub>&sup2; ) )

  where the symmetry axis is the Cu-Cu direction for Cu and the O&rarr;Cu bond direction for O.
- A fraction &rho;<sub>Cu</sub> of the spin sits on each Cu and (1 - &rho;<sub>Cu</sub>)/4 on each of its four oxygens.
- The cluster structure factor is S<sub>c</sub>(Q) = &sum; amp &middot; f(Q) &middot; exp(i Q&middot;r). For a dimer (a, b) the triplet-excitation intensity is proportional to |S<sub>a</sub> - S<sub>b</sub>|&sup2;, and the two dimers add incoherently:

  I(Q) = scale &middot; ( |S<sub>1</sub> - S<sub>2</sub>|&sup2; + |S<sub>3</sub> - S<sub>4</sub>|&sup2; ) + background

  No polarisation factor or separate magnetic form factor is applied.
- Because the data are voxel averages, the model is averaged over `nsub`&sup3; points per voxel.

### Free parameters

| # | Name | Meaning | Unit |
|---|------|---------|------|
| 1 | `rho_Cu` | spin fraction on each Cu (0 to 1) | - |
| 2 | `sigma_Cu_perp` | Cu Gaussian width perpendicular to the Cu-Cu axis (within the CuO<sub>4</sub> plane) | &Aring; |
| 3 | `sigma_Cu_par` | Cu Gaussian width along the Cu-Cu axis | &Aring; |
| 4 | `sigma_O_perp` | O Gaussian width perpendicular to the O-Cu bond | &Aring; |
| 5 | `sigma_O_par` | O Gaussian width along the O-Cu bond | &Aring; |
| 6 | `scale` | overall intensity scale | arb. |

The flat background is fixed (default 0, `cfg.background`).

## Requirements

- MATLAB R2022a or newer (R2017b is enough if you drop the plots that use `clim`)
- Optimization Toolbox (`lsqnonlin`)
- Parallel Computing Toolbox (optional; `parfor` runs serially without it)

## Data

Raw data are available from https://doi.ill.fr/10.5291/ILL-DATA.4-05-879, emabrgo period my apply.

Place the following `.mat` files in a `data/` folder next to the script. `<dQ>` is the voxel size, 0.15 by default.

| File | Variable | Content |
|------|----------|---------|
| `I_background_subtracted_dQ_<dQ>.mat` | `sub` | 3D intensity, ordered (Qx, Qy, Qz) |
| `I_err_background_subtracted_dQ_<dQ>.mat` | `sub_err` | its uncertainty |
| `Qx_background_subtracted_dQ_<dQ>.mat` | `Qx` | grid axis, used as H |
| `Qy_background_subtracted_dQ_<dQ>.mat` | `Qy` | grid axis, used as L |
| `Qz_background_subtracted_dQ_<dQ>.mat` | `Qz` | grid axis, used as K |

Voxels with zero intensity or zero error are treated as empty and excluded from the fit.

### Coordinates

The grid axes (H, K, L) are **not** Miller indices. A grid point corresponds to the reciprocal-lattice vector

```
hkl = H*u1 + L*u2 + K*[0 1 0]      (r.l.u.)
```

where `u1` and `u2` (set in the configuration) are the in-plane directions of the data grid and K runs along b\*. The Cartesian momentum transfer in &Aring;<sup>-1</sup> is `Q = B*hkl'` with `B = 2*pi*inv(A)`, where the rows of `A` are the lattice vectors (monoclinic cell, a along x, b along y). Plot axes are labelled with the true hkl direction and the length of one unit in &Aring;<sup>-1</sup>.

## Usage

1. Put the data in `data/` (or change `cfg.dataDir`).
2. Open `fit_spin_density_Cu_acetate.m` and check the configuration section at the top.
3. Run the script.

The fit prints a table of parameters and uncertainties. Uncertainties come from the Jacobian, rescaled by the reduced chi-squared, and are not meaningful for a parameter that ends on a bound. If `cfg.resultsFile` is set (default `fit_results.mat`), the results are saved.

### Main settings

| Setting | Meaning |
|---------|---------|
| `cfg.dQ` | voxel size of the data grid (also part of the file names) |
| `cfg.u1`, `cfg.u2` | in-plane directions of the grid in r.l.u. |
| `cfg.crop` | region of the grid that is loaded |
| `cfg.mask` | regions excluded from the fit |
| `cfg.nsub` | voxel sub-sampling, `nsub`&sup3; points per voxel |
| `cfg.fit.p0`, `lb`, `ub` | start values and bounds |
| `cfg.cut` | positions of the 1D cuts |



## Repository contents

```
fit_spin_density_Cu_acetate.m   main script (configuration, fit, plots, local functions)
CITATION.cff                    citation metadata
LICENSE                         licence
```

## Citation

If you use this code, please cite the article (see `CITATION.cff` or the GitHub "Cite this repository" button):

Federico Pisani, Leonie Spitz, Libor Vojáček et al., Probing the real-space density of spin-entangled electrons, 2026, arXiv:2604.13751.

## Licence

Released under the MIT licence, see `LICENSE`.

## Contact

Federico Pisani
Laboratory for quantum magnetism, EPFL
federico.pisani@epfl.ch
