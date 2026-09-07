# Charge-density difference (Δρ) — hero drug/2D-material complex

`Δρ(r) = ρ_complex(r) − ρ_carrier(r) − ρ_drug(r)`

All three electron densities are evaluated with **GFN2-xTB** at the *complex*
geometry (the drug and carrier fragments keep their bound coordinates), sampled
on one shared grid. Yellow isosurface lobes = electron accumulation, blue =
depletion. This visualises the interfacial charge reorganisation that
accompanies the adsorption energy reported in the main text.

## Files
| file | description |
|---|---|
| `*_complex.xyz` | GFN2-xTB optimised drug/2D-material complex (from `calculations/`) |
| `*_carrier.xyz` | 2D-material fragment at complex coordinates |
| `*_drug.xyz` | drug fragment at complex coordinates |
| `*_deltarho.cub.gz` | gzipped Gaussian cube of Δρ |
| `*_deltarho_render.png` | ChimeraX render used in the figure |
| `build_deltarho.py` | regenerates the cube from the three `.xyz` (needs xtb + Multiwfn 3.8) |

## Reproduce
```bash
cd results/quantum/drho
python build_deltarho.py <key>      # <key> = kras | gbm | tau | tnbc
```
Fragments were separated by chemical identity (2D-material marker element vs the
rest); each H atom assigned to its nearer heavy atom. No geometry was altered.

The figure is assembled by `src/visualization/_drho_fig.py`, called from the
master figure generator.
