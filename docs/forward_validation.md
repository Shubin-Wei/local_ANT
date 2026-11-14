# The 1-D forward calculation (surfdisp96)

[中文](forward_validation.zh.md) | English

Every phase velocity of a layered model in lant is computed by DSurfTomo's own dispersion code: the starting model and its fit, the kernel figure, the synthetic examples and, column by column, the depth kernels of the 3-D inversion. `lant.inversion.dsurf_disp96` ports DSurfTomo `a651faf`'s `surfdisp96.f` (Rayleigh phase velocity, spherical earth, fundamental mode) and, from `CalSurfG.f90`, the Brocher materials and `refineGrid2LayerMdl`, statement by statement and in the Fortran's own precisions. The phase velocities and the depth kernels equal those of the native `depthkernel` bit for bit (a stored native output is part of the tests), and the 3-D inversion is byte-identical to native DSurfTomo on four test cases ([python_backends.md](python_backends.md)).

Before 2.0 the 1-D calculation was first fwd1d, a small Fortran program around surfdisp96, then disba 0.7.0 (surf96 in double precision); disba's phase velocities differed from surfdisp96's by up to 5e-6 km/s, enough to move single data across DSurfTomo's rejection bound.

## API and model semantics

```python
from lant.inversion.start import forward

c = forward([0, 0.1, 0.5, 2], [1.8, 2.0, 2.6, 3.4], [1, 0.1, 0.5])
```

Depth in km, Vs and the result in km/s, periods in s. Only Rayleigh phase velocities (`igr=0`, `iwave` 0 or 2) are ported. At least two depth nodes starting at zero and increasing, positive Vs, and 1–80 positive periods (NP of surfdisp96) are required. The periods are solved in increasing order, as DSurfTomo's are (surfdisp96 starts the root search of each period from the previous root), and returned in the caller's order. A period without a root is returned as 0, as surfdisp96 does; the callers stop on it.

As in DSurfTomo, the node depths and Vs are REAL (single precision). Vp and density follow Brocher (2005) at the nodes and are interpolated to the sublayer midpoints. Each interval of thickness `h` becomes `int((h + 1e-4) / (h / minthk)) + 1` sublayers; with the half-space at most 200 layers. Spherical flattening uses radius 6370 km, logarithmic thickness, midpoint velocity scaling and the density power −2.275 (Rayleigh), with the half-space scaled at a temporary 1 km thickness. The root search is surfdisp96's: a starting velocity of 0.95 × 0.90 of the half-space Rayleigh speed of the slowest layer, steps of `dc = 0.005` (as a REAL constant: 0.004999999888 km/s), Dunkin's compound-matrix dispersion function in double precision, refinement by Neville interpolation and bisection to 1e-6 relative; the result is stored as REAL.

## Limitations

Higher modes, group velocities, Love waves, fluid layers and anisotropy are not ported; lant uses fundamental Rayleigh phase velocities only. Where the root search fails, native DSurfTomo prints "improper initial value in disper" and continues with 0; lant stops with that message.

## Sources

- [DSurfTomo a651faf](https://github.com/HongjianFang/DSurfTomo/tree/a651faf92758c1b4c5eb494a4865ef797b9a4b1c), MIT: `surfdisp96.f` and `CalSurfG.f90` (`depthkernel`, `caldespersion`, `refineGrid2LayerMdl`).
- Herrmann (2013), *Seismol. Res. Lett.* 84(6), 1081–1088, doi:10.1785/0220110096: Computer Programs in Seismology, whose surf96 `surfdisp96.f` modifies.
- Brocher (2005), *BSSA*, [10.1785/0120050077](https://doi.org/10.1785/0120050077): the empirical Vp/Vs/density relationships of DSurfTomo.
- Schwab & Knopoff (1972), *Methods in Computational Physics*, vol. 11, and Biswas (1972), *PAGEOPH* 96, 61–74: the spherical flattening and density mapping, as cited in surfdisp96's `sphere`.
