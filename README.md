# Competition Rifle Ballistics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia Version](https://img.shields.io/badge/julia-1.6+-blue.svg)](https://julialang.org)

A comprehensive mathematical library and toolset for long-range competition rifle ballistics, implemented in Julia. This project provides everything from basic unit conversions to high-fidelity 6-Degree-of-Freedom (6-DOF) trajectory simulations.

This repository is the companion code for the manual: ***"Competition Rifle Ballistics: A Comprehensive Mathematical Manual"***.

## 🚀 Overview

`CompetitionBallistics.jl` is designed for competitive shooters, handloaders, and ballistics researchers who require high-precision trajectory modeling and data-driven reloading analysis.

### Key Features

*   **3-DOF Trajectory Solver:** High-precision RK4 integrator accounting for gravity, drag, Coriolis effect (horizontal and vertical/Eötvös), and spin drift.
*   **6-DOF Rigid-Body Dynamics:** Advanced solver for angular motion (yaw, precession, nutation) and first-principles spin drift prediction.
*   **Atmospheric Modeling:** Full ICAO standard atmosphere with humidity corrections and density altitude calculations.
*   **Drag Models:** Support for G1/G7 standard tables with C²-continuous cubic spline interpolation, plus support for custom Doppler-radar measured drag profiles.
*   **Interior Ballistics:** Le Duc velocity/pressure models, barrel length corrections, and powder temperature sensitivity analysis.
*   **Reloading Analysis:** Tools for Satterlee ladder tests, chronograph statistics (ES/SD), vertical dispersion estimation, and BC estimation from dual-chronograph data.
*   **Interactive Web App:** A built-in Dash.jl dashboard for visual trajectory analysis.
*   **Pluto Tutorial:** A step-by-step interactive notebook for learning ballistic principles.

## 📂 Project Structure

```text
.
├── src/
│   ├── CompetitionBallistics.jl  # Core library module (Julia — the reference)
│   ├── CompetitionBallistics.js  # JavaScript port — the copy served by tireur.org
│   ├── app.jl                   # Dash.jl web application
│   └── BallisticsTutorial.jl     # Pluto.jl tutorial notebook
├── doc/
│   └── ballistics_manual.tex     # LaTeX source for the manual
├── LICENSE                      # MIT License
└── README.md                    # You are here
```

The two ports are kept in line by a nightly validation harness living in the
**parent repository** ([tireur.org](https://github.com/fbastin/tireur.org), where
this repo is vendored as a submodule) — see
[Replaying the validation controls](#replaying-the-validation-controls).

## 🛠️ Installation & Setup

### Prerequisites

1.  **Julia:** Install Julia (v1.6 or later) from [julialang.org](https://julialang.org/downloads/).
2.  **Dependencies:** Open the Julia REPL and install the required packages:

```julia
using Pkg
Pkg.add(["Dash", "PlotlyJS", "Pluto", "Statistics", "Printf", "LinearAlgebra"])
```

## 📖 Usage

### Using the Library

To use the ballistics engine in your own Julia scripts:

```julia
include("src/CompetitionBallistics.jl")
using .CompetitionBallistics.ExteriorBallistics

# Define a shot (6.5 Creedmoor example)
params = ShotParameters(
    mass_grains = 140.0,
    caliber_in = 0.264,
    bc = 0.311,
    drag_model = :G7,
    muzzle_vel_fps = 2700.0,
    target_range_yd = 1000.0
)

# Solve and print a DOPE table
trajectory_table(params, step_yd=100)
```

### Running the Web App

Launch the interactive dashboard to visualize trajectories and compare bullets:

```bash
julia src/app.jl
```
Then open your browser at `http://localhost:8050`.

### Interactive Tutorial

For a guided walkthrough of the library's features, launch the Pluto notebook:

```julia
using Pluto
Pluto.run(notebook="src/BallisticsTutorial.jl")
```

## ✅ Replaying the Validation Controls

Three independent controls watch this library. None trusts the others — they are
built to catch *different* failure modes:

1. **Parity (Julia ↔ JavaScript)** — the two ports are separate hand-written
   implementations of the same manual. They are run on the same grid of inputs
   (86 cases over 21 functions, plus 45 reference data values) and compared at
   1e-9 relative tolerance. Catches port drift. `--selftest` perturbs one
   ballistic coefficient by 1 % and verifies the guard bites.
2. **Frozen third-party reference (JavaScript ↔ `py-ballisticcalc`)** — parity
   cannot see an error *shared by both ports* (this happened: identical
   supersonic drag tables, wrong above Mach 1.03, with perfect parity). The
   solver is therefore confronted with a frozen, provenance-stamped reference
   produced by the independent library `py-ballisticcalc` (point-masse 3-DOF,
   RK4 engine): 14 cases, 162 points. A frozen file, not a live call — a
   reference hosted elsewhere is a revocable loan, not a fixture.
3. **Analytic cases and convergence** — standard atmosphere density, speed of
   sound, barometric pressure, zero-angle search, form factor identities.

### Prerequisites

* **Julia** ≥ 1.6 (no extra packages: the controls only use the standard library)
* **Node.js** (any recent version)
* **Python 3** (standard library only)
* To *regenerate* the frozen reference (optional, deliberate): `pip install py-ballisticcalc`

### Running

```bash
# clone the parent repository — it carries this repo as a submodule and the harness
git clone --recurse-submodules https://github.com/fbastin/tireur.org
cd tireur.org

# 1. parity between the two ports (fails on any drift)
python3 scripts/check_js_parity.py            # add --selftest to verify the guard bites

# 2. our JavaScript solver vs the frozen third-party reference
python3 scripts/check_ballistics_reference.py --selftest   # guard self-check
python3 scripts/check_ballistics_reference.py              # 14 cases, 162 points

# 3. (optional, deliberate) regenerate the frozen reference — it is provenance-stamped
#    (library version, engine, date) and must not be refreshed casually:
python3 scripts/gen_ballistics_reference.py   # requires py-ballisticcalc
```

All three print `OK` on a healthy checkout. If a divergence appears, the rule is:
**arbitrate against the published tables (McCoy, Litz) or measured data — never
align one implementation to the other.** The frozen reference is a second opinion,
not an oracle; two implementations can share the same error (ours did, for months).

## 📚 Mathematical Manual

The theoretical foundations, including all derivations for the 3-DOF and 6-DOF models, are detailed in the `doc/ballistics_manual.tex` file. You can compile this with `pdflatex` to generate the companion PDF.

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
*Created by Fabian Bastin (2026)*
