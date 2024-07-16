# Tools for producing the Muon Collider benchmark plots

Collection of configuration files and scripts for obtaining the basic overview of the Muon Collider detector performance, also used for software-release validation.

## Repository structure

All the code is split into separate stages of the performance-evaluation process:

* `generation` - production of a sample of input objects, typically `MCParticles`;
* `simulation` - simulation of particles' interaction with the detector in GEANT4, typically producing `SimHits`;
* `digitisation` - application of detector-level effects and filtering to the `SimHits`, typically producing `RecHits`;
* `reconstruction` - reconstruction of high-level objects, such as `Tracks`, `Jets`, `ParticleFlow` objects, etc.;
* `analysis` - production of files with standardised objects types compatible with generic plotting scripts, typically `ROOT.TTree`, `ROOT.TH1`, etc.
* `plotting` - generic scripts for producing plots from the `analysis` output.

Specific workflow scripts and corresponding plotting configurations for individual use cases are stored under `workflows/` directory.
Each workflow can use any subset of the `generation` to `plotting` stages, depending on its purpose.

