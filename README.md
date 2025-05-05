# SwarmMTFilGliding

> **Gliding assay at high densities of microtubules**  
> In certain in vitro gliding assays (e.g., Sumino et al.), dense layers of microtubules exhibit spontaneous vortex, milling and flocking patterns.  
> This repository provides minimal Julia implementations of four active‐matter models designed to capture and explore these collective phenomena.

---

## Repository Structure

```text
.
├── Blind-Angle/
│   └── Blind-Angle.jl        # Vicsek-like model with a “blind field of view” behind each particle
│
├── Nematic-Rods/
│   └── nematic-rods.jl        # Self‐propelled rods with nematic (headless) alignment
│
├── SPP-Vicsek/
│   └── vicsek-classic.jl      # Classic Vicsek flocking model (polar alignment)
│
└── SPP-with-memory/
    └── spp-memory.jl          # Memory‐enhanced SPP model with intrinsic angular inertia
```

Each folder contains:
- A single `*.jl` script implementing the model.
- A single `*.ipynb` notebook for analysis of the generated CSV files.

---

## Getting Started

### Prerequisites

- **Julia 1.8+**  
- Standard libraries:
  ```julia
  using Random, Statistics, Plots
  ```
- (Optionally) increase threading:
  ```bash
  export JULIA_NUM_THREADS=4
  ```

### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/CyCelsLab/swarmMTFilGliding.git
   cd swarmMTFilGliding
   ```
### Running a Simulation

From the repository root, navigate to the model directory and run the script. For example:

```bash
cd Blind-Angle
julia Blind-Angle.jl
```

This produces:
- **Particle data** CSV: `particles_L{L}_… .csv`
- **Order parameters** CSV: `order_L{L}_… .csv`
- **(Memory model only)** PNG snapshots: `particles_t{t}.png`

Parameters (domain size L, density ρ, noise η, etc.) can be edited at the top of each script before running.

---

## Model Summaries

| Model            | Alignment Type        | Memory/Noise      | Key Reference                                       |
| ---------------- | --------------------- | ----------------- | ----------------------------------------------------|
| **Vicsek**       | Polar                 | Angular noise     | Vicsek *et al.* (1995) `Phys. Rev. Lett. 75, 1226`   |
| **Blind‐Angle**  | Polar + blind FoV     | Angular noise     | Costanzo & Hemelrijk (2018) `J. Phys. D: Appl. Phys. 51` |
| **Nematic‐Rods** | Nematic               | Angular noise     | Ginelli *et al.* (2010) `Phys. Rev. Lett. 104, 184502` |
| **Memory SPP**   | Nematic + Memory      | Noise + τ memory  | Sumino *et al.* (2012) `Nature 483, 448` <br>Nagai *et al.* (2015) `Phys. Rev. Lett. 114, 168001` |
