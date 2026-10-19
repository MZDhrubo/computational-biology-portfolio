# Molecular-dynamics reproducibility checklist

My CV records GROMACS-based molecular-dynamics training. This checklist is a learning and reporting aid, not evidence that a complete production simulation was performed as part of this repository.

For a future or separately documented simulation, record:

- GROMACS version and operating environment;
- provenance and preparation of the starting structure;
- force field, water model, ions, box dimensions, and topology-generation choices;
- energy-minimization and equilibration settings;
- ensemble, temperature/pressure coupling, timestep, constraints, and production length;
- random seeds and configuration files where available;
- all warnings, convergence checks, and failed/restarted runs;
- trajectory-processing selections and analysis scripts;
- limits of the simulation and the distinction between computational observations and experimental evidence.

Never publish private structures, collaborator data, unpublished trajectories, or project-specific results without permission. A reproducible setup is necessary but does not by itself establish that the simulated system is biologically accurate.
