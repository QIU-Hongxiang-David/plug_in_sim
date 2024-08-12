# plug_in_sim

Code to replicate the simulation in "Non-Plug-In Estimators Could Outperform Plug-In Estimators: a Cautionary Note and a Diagnosis".

- `main_text/`: Simulation in the main text and for the TMLE with a small known bound on the outcome regression.
- `large_Q_bound/`: Simulation for the TMLE with a relatively large known bound on the outcome regression.

Under each folder:

- `setup.R`: simulation setup. The `run.once` function contains code to run one simulation repetition.
- `simulation.R`: R code to run the simulation. Designed to run on a computing cluster with R package `batchtools` set up properly.
- `results.R` (only under `main_text/`): R code to generate plots and tables for both simulations.
