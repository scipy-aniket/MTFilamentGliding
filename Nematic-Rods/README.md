# Nematic Rods Simulation in Julia

This repository contains a Julia implementation of the large-scale collective behavior of self-propelled rods with nematic alignment, based on the work of Ginelli *et al.* 2010. The simulation captures emergent nematic ordering as a function of density, noise, and interaction radius.

## File Structure

* `nematic-rods.jl`: Main simulation script.
* Output files:

  * `particles_L{L}_T{T}_N{N}_eta{η}.csv`: Particle positions and orientations over time.
  * `order_parameters_L{L}_T{T}_N{N}_eta{η}.csv`: Polar and nematic order parameters over time.

## Overview

In this model, point-like rods move with constant speed in a two-dimensional periodic domain. Each rod aligns its orientation nematically (headless) with neighbors within an interaction radius, perturbed by angular noise. Spatial partitioning into grid cells accelerates neighbor searches.

## Parameters

Defined at the top of `nematic-rods.jl`:

| Parameter            | Description                                    | Example     |
| -------------------- | ---------------------------------------------- | ----------- |
| `L`                  | Domain side length                             | 30          |
| `T`                  | Number of timesteps                            | 1000        |
| `ρ` (rho)            | Particle density                               | 0.5         |
| `v_0`                | Speed of each rod                              | 0.5         |
| `cell_length`        | Grid cell size                                 | 1           |
| `interaction_radius` | Neighborhood radius for nematic alignment      | 1.0         |
| `η` (eta)            | Noise amplitude                                | 0.1         |
| `noise_lims`         | White noise limits (± range for angular noise) | (-π/2, π/2) |
| `dt`                 | Time step size                                 | 1.0         |
| `frameskip`          | Output data every this many steps              | 2           |
| `N`                  | Number of particles (computed as `Int(ρ*L^2)`) | —           |

## Running the Simulation

Ensure Julia is installed, then from the project directory:

```bash
julia nematic-rods.jl
```

This generates CSV files in the current directory.

## Dependencies

The script requires the following Julia packages:

* `Random`
* `Plots`
* `Statistics`

Install them via:

```julia
using Pkg
Pkg.add(["Plots", "Statistics"])
```

## Output

1. **Particle Data (`particles_*.csv`)**:

   * Columns: `time_step`, `particle_id`, `x`, `y`, `theta`
   * State of each rod at each output frame

2. **Order Parameters (`order_parameters_*.csv`)**:

   * Columns: `time_step`, `P_t`, `S_t`
   * `P_t`: Polar order parameter (with head orientation)
   * `S_t`: Nematic order parameter (headless alignment)

## Performance Notes

* Uses multi-threading (`Threads.@threads`) for updating rod states.
* Spatial grid lookup reduces neighbor search complexity.
* Control threads via:

```bash
JULIA_NUM_THREADS=4 julia nematic-rods.jl
```

## Citation

Based on:

Ginelli, Francesco, Fernando Peruani, Markus Bär, and Hugues Chaté. *"Large-scale collective properties of self-propelled rods."* Physical Review Letters 104, no. 18 (2010): 184502.

BibTeX:

```bibtex
@article{ginelli2010large,
  title={Large-scale collective properties of self-propelled rods},
  author={Ginelli, Francesco and Peruani, Fernando and B{"a}r, Markus and Chat{'e}, Hugues},
  journal={Physical review letters},
  volume={104},
  number={18},
  pages={184502},
  year={2010},
  publisher={APS}
}
```
