# Purpose-specific workflows

Each workflow comprises of two components:
* the main workflow script `workflow.sh` that calls sequentially all the relevant stages: `generation`, ..., `plotting`;
* supporting configuration files for the generic `analysis` and `plotting` scripts tailored to the specific workflow.

## Available workflows

* `relval` - run a full cycle from `generation` to `reconstruction` on limited statistics for software-release validation;
* `bib_production` - run GEANT4 simulation of a complete BIB sample as part of central sample production;
* `mdi_analysis` - run GEANT4 simulation and reconstruction of a small BIB sample as part of MDI optimisation;
