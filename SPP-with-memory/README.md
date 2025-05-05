# Memory-Enhanced SPP Simulation in Julia

This repository contains a Julia implementation of a self-propelled particle (SPP) model with memory effects, based on the work of Sumino *et al.* (2012) and Nagai *et al.* (2015). Particles carry an internal angular velocity (ω) that evolves with memory and noise, leading to complex collective patterns such as vortex lattices and swirling.

## File Structure

* `spp-memory.jl`: Main simulation script.
* Output files:

  * `L{L}_N{N}_T{T}_ρ{ρ}_w0{ω_0}_s{σ_noise}_a{α}_t{τ}.csv`: Time series of positions, orientations, and angular velocities.
  * `particles_t{t}.png`: Snapshot plots of particle positions every 200 steps.

## Overview

In this model, N particles move at constant speed in a periodic 2D box. Each particle carries an intrinsic angular velocity that relaxes toward ω₀ with memory time τ and is subject to Gaussian noise of strength σ. Nematic interactions with neighbors within a fixed radius induce torques, coupling orientation updates to both neighbor alignment and memory.

## Parameters

Defined at the top of `spp-memory.jl`:

| Parameter            | Description                          | Example |
| -------------------- | ------------------------------------ | ------- |
| `L`                  | Domain side length                   | 64      |
| `T`                  | Number of timesteps                  | 1200    |
| `ρ` (rho)            | Particle density                     | 5       |
| `v_0`                | Particle speed                       | 0.5     |
| `cell_length`        | Grid cell size for neighbor lists    | 1       |
| `interaction_radius` | Neighborhood radius for interactions | 1       |
| `dt`                 | Time step size                       | 1       |
| `N`                  | Number of particles (`Int(ρ*L^2)`)   | —       |
| `ω_0`                | Basal angular velocity               | −0.1    |
| `σ_noise`            | Noise standard deviation on ω update | 0.003   |
| `α`                  | Strength of nematic torque           | 0.1     |
| `τ` (tau)            | Memory time scale for ω relaxation   | 100     |
| `frameskip`          | Output interval for CSV data         | 5       |

## Running the Simulation

Ensure Julia is installed, then run:

```bash
julia spp-memory.jl
```

This produces:

* A CSV file containing `time_step,particle_id,x,y,theta,omega` for each frameskip.
* PNG snapshots `particles_t0.png`, `particles_t200.png`, etc., showing spatial configurations.

## Dependencies

The script requires:

* `Random`
* `Plots`
* `Statistics`

Install with:

```julia
using Pkg
Pkg.add(["Plots", "Statistics"])
```

## Output

1. **Time-Series Data** (`.csv`):

   * Columns: `time_step,particle_id,x,y,theta,omega`
   * Records positions, orientations, and angular velocities at each recorded step.

2. **Snapshots** (`.png`):

   * Particle position scatter plots every 200 time steps.
   * File names: `particles_t{t}.png`.

## Performance Notes

* Spatial grid lookup accelerates neighbor searches.
* Multi-threaded pre-update and neighbor update via `Threads.@threads`.
* Adjust thread count:

```bash
JULIA_NUM_THREADS=4 julia spp-memory.jl
```

## Citation

Sumino, Y., Nagai, K. H., Shitaka, Y., Tanaka, D., Yoshikawa, K., Chaté, H., & Oiwa, K. *"Large-scale vortex lattice emerging from collectively moving microtubules."* Nature 483, 448–452 (2012).

Nagai, K. H., Sumino, Y., Montagne, R., Aranson, I. S., & Chaté, H. *"Collective motion of self-propelled particles with memory."* Phys. Rev. Lett. 114, 168001 (2015).

BibTeX:

```bibtex
@article{sumino2012large,
  title={Large-scale vortex lattice emerging from collectively moving microtubules},
  author={Sumino, Yutaka and Nagai, Ken H. and Shitaka, Yuji and Tanaka, Dan and Yoshikawa, Kenichi and Chat{\'e}, Hugues and Oiwa, Kazuhiro},
  journal={Nature},
  volume={483},
  number={7390},
  pages={448--452},
  year={2012},
  publisher={Nature Publishing Group}
}

@article{nagai2015collective,
  title={Collective motion of self-propelled particles with memory},
  author={Nagai, Ken H. and Sumino, Yutaka and Montagne, Raul and Aranson, Igor S. and Chat{\'e}, Hugues},
  journal={Physical Review Letters},
  volume={114},
  number={16},
  pages={168001},
  year={2015},
  publisher={APS}
}
```
