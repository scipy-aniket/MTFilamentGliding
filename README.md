# Blind-Angle Simulation in Julia

This repository contains a Julia implementation of a Vicsek-like model incorporating a blind spot (blind angle) for each self-propelled particle, based on the work of Costanzo & Hemelrijk (2018). The simulation explores how a blind angle affects collective patterns such as milling and flocking.

## File Structure

* `Blind-Angle.jl`: Main simulation script.
* Output files:

  * `modified_vicsek_particles_L{L}_T{T}_N{N}_eta{η}_omega{ω}.csv`: Particle positions and orientations over time.
  * `modified_vicsek_order_L{L}_T{T}_N{N}_eta{η}_omega{ω}.csv`: Polar, nematic order parameters, and angular momentum over time.

## Overview

In this model, particles move at constant speed in a 2D periodic box. Each particle aligns its direction with visible neighbors within an interaction radius, excluding those within a specified blind angle behind it. Angular noise and a maximum turning rate (ω) further modulate the dynamics.

## Parameters

Defined at the top of `Blind-Angle.jl`:

| Parameter            | Description                                      | Example      |
| -------------------- | ------------------------------------------------ | ------------ |
| `L`                  | Domain side length                               | 20           |
| `T`                  | Number of timesteps                              | 1000         |
| `ρ` (rho)            | Particle density                                 | 2.0          |
| `v_0`                | Particle speed                                   | 0.1          |
| `cell_length`        | Grid cell size for neighbor lists                | 1.0          |
| `interaction_radius` | Neighborhood radius for interaction              | 1.0          |
| `dt`                 | Time step size                                   | 1.0          |
| `N`                  | Number of particles (`Int(ρ*L^2)`)               | —            |
| `ω` (omega)          | Maximal Turning Rate (`v_0/(r*2.03)`)            | computed     |
| `η` (eta)            | Noise strength (`0.1*ω`)                         | computed     |
| `noise_bounds`       | Bounds for angular noise                         | (-0.5, 0.5)  |
| `blind_angle`        | Range of angles defining blind spot (in radians) | \[10°, 350°] |
| `frameskip`          | Output data every this many steps                | 2            |

## Running the Simulation

Ensure Julia is installed. From the project directory:

```bash
julia Blind-Angle.jl
```

This will generate two CSV files:

* `modified_vicsek_particles_...csv`
* `modified_vicsek_order_...csv`

## Dependencies

The script uses:

* `Random`
* `Statistics`

Install via Julia's package manager if needed:

```julia
using Pkg
Pkg.add(["Statistics"])
```

## Output

1. **Particle Data**:

   * Columns: `time_step,particle_id,x,y,theta`
   * Positions and orientations at each output frame.

2. **Order Parameters**:

   * Columns: `time_step,P_t,S_t,m_a`
   * `P_t`: Polar order.
   * `S_t`: Nematic order.
   * `m_a`: Mean absolute angular momentum about center of mass.

## Performance Notes

* Uses spatial grid lookup to speed neighbor searches.
* Multi-threaded updates through `Threads.@threads`.
* Control threads with:

```bash
JULIA_NUM_THREADS=4 julia Blind-Angle.jl
```

## Citation

Based on:

Costanzo, A., & Hemelrijk, C.K. *"Spontaneous emergence of milling (vortex state) in a Vicsek-like model."* Journal of Physics D: Applied Physics 51, no. 13 (2018): 134004.

BibTeX:

```bibtex
@article{costanzo2018spontaneous,
  title={Spontaneous emergence of milling (vortex state) in a Vicsek-like model},
  author={Costanzo, A and Hemelrijk, CK},
  journal={Journal of Physics D: Applied Physics},
  volume={51},
  number={13},
  pages={134004},
  year={2018},
  publisher={IOP Publishing}
}
```
