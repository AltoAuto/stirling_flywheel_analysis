# Stirling Flywheel Analysis

MATLAB implementation of a **beta-type Stirling engine analysis and flywheel design**, developed for **ME 4053: Energy Conversion & Storage (Fall 2025, University of Minnesota)**.

The project models a Stirling engine with crank–slider kinematics for both the displacer and power piston, computes torque and power output, and sizes a flywheel that meets the design specifications in Appendix B. The analysis is fully parametric and produces all required plots automatically.

---

## Features

* **Single-file MATLAB script (`main.m`)**: all helpers are local functions for easy grading.
* **Configurable assumptions**: ideal gas, isothermal compression/expansion, perfect regenerator.
* **End-to-end workflow**:

  * Kinematics → volumes vs crank angle
  * Thermodynamic pressure model (Schmidt-style)
  * Torque & power (two independent methods)
  * Energy deviation & required flywheel inertia
  * Speed ripple simulation
  * Phase-angle sweep
* **Auto-generated outputs**:

  * Pressure vs. specific volume (engine vs. ideal cycle)
  * Torque vs. crank angle
  * Rotational velocity vs. crank angle
  * Energy per cycle vs. phase angle
  * Summary table with power, inertia, and flywheel dimensions

---

## Quickstart

Clone/download this repo and run the script in MATLAB:

```matlab
clear; close all; clc;
main
```

The script will:

1. Print configuration, assumptions, and summary results to the console.
2. Generate required plots and save them in `outputs/`.
3. Write a `summary.csv` file with key numbers (power, inertia, flywheel specs).

---

## Repo Structure

```
stirling-flywheel-analysis/
│
├── main.m          # one-file MATLAB implementation
├── outputs/        # auto-saved plots + summary.csv
├── report/         # executive summary draft & figures
└── README.md       # project overview (this file)
```

---

## Deliverables (per project brief)

* Basic flywheel design and dimensions
* Power output (via two methods)
* Required plots:

  * PV diagram (Stirling cycle + engine loop)
  * Torque vs. crank angle
  * Rotational velocity vs. crank angle
  * Energy per cycle vs. phase angle
* Textual analysis and summary (report, 6–10 pages)

---

## References

* Project description & specifications:

  * *Stirling Engine Analysis & Flywheel Design*
  * Appendix A: Project Expectations
  * Appendix B: Problem Specifications
