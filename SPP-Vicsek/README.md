
# Vicsek Model Simulation in Julia

This repository contains a Julia implementation of the classic Vicsek model, which simulates the collective motion of self-propelled particles (SPPs) interacting locally in two dimensions. The model demonstrates emergent behavior such as flocking and phase transitions as a function of noise and density.

## File Structure

- `SPP-Vicsek/vicsek-classic.jl`: Main simulation script.
- Output files:
  - `particles_L{L}_T{T}_N{N}_eta{η}.csv`: Particle positions and orientations over time.
  - `order_parameters_L{L}_T{T}_N{N}_eta{η}.csv`: Polar and nematic order parameters over time.

## Overview

The Vicsek model consists of point particles moving with constant speed. Each particle aligns its direction with the average orientation of its neighbors within a fixed radius, perturbed by some angular noise. The model uses periodic boundary conditions and spatial partitioning into grid cells to accelerate neighbor search.

## Parameters

Defined at the top of the script:

| Parameter           | Description                                      | Example |
|---------------------|--------------------------------------------------|---------|
| `L`                 | Length of simulation domain                      | 30      |
| `T`                 | Number of timesteps                              | 1000    |
| `ρ` (rho)           | Particle density                                 | 0.5     |
| `v_0`               | Speed of each particle                           | 0.5     |
| `interaction_radius`| Neighborhood radius for alignment                | 1.0     |
| `η` (eta)           | Noise amplitude                                  | 0.5     |
| `dt`                | Time step size                                   | 1.0     |
| `frameskip`         | Output data every `frameskip` steps              | 2       |

The number of particles `N` is computed as `N = ρ * L^2`.

## Running the Simulation

To run the simulation, make sure Julia is installed, then from the project directory run:

```bash
julia SPP-Vicsek/vicsek-classic.jl
```

This will generate the CSV output files in the current working directory.

## Dependencies

The simulation script uses the following Julia packages:

- `Random`
- `Plots`
- `Statistics`

You can install the required packages using:

```julia
using Pkg
Pkg.add(["Plots", "Statistics"])
```

## Output

Two CSV files are generated:

1. **Particle Data (`particles_*.csv`)**:
    - Columns: `time_step`, `particle_id`, `x`, `y`, `theta`
    - Contains the state of each particle at each output frame

2. **Order Parameters (`order_parameters_*.csv`)**:
    - Columns: `time_step`, `P_t`, `S_t`
    - `P_t`: Polar order parameter (mean alignment)
    - `S_t`: Nematic order parameter (headless alignment)

## Performance Notes

- The simulation uses multi-threading via `Threads.@threads` for updating particles.
- You can control the number of threads by setting the environment variable `JULIA_NUM_THREADS` before running:

```bash
JULIA_NUM_THREADS=4 julia SPP-Vicsek/vicsek-classic.jl
```

## Citation

This implementation is based on the original Vicsek model introduced in:

Vicsek, Tamás, et al. *"Novel type of phase transition in a system of self-driven particles."* Physical review letters 75.6 (1995): 1226.

BibTeX:
```bibtex
@article{vicsek1995novel,
  title={Novel type of phase transition in a system of self-driven particles},
  author={Vicsek, Tam{\'a}s and Czir{\'o}k, Andr{\'a}s and Ben-Jacob, Eshel and Cohen, Inon and Shochet, Ofer},
  journal={Physical review letters},
  volume={75},
  number={6},
  pages={1226},
  year={1995},
  publisher={APS}
}
```

## License

This code is provided for academic and research purposes. Modify and use it at your own discretion.
