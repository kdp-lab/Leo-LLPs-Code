# Leo-LLPs-Code

This repository contains code for running multi-processing over sim, digi, and reco for studying long-lived particles (LLPs).


### Prerequisites:
- Clone both [`MuC-tutorial`](https://github.com/MuonColliderSoft/MuC-Tutorial) and [`mucoll-benchmarks`](https://github.com/MuonColliderSoft/mucoll-benchmarks/tree/main) into this directory

## Instructions for running multi_sim.py:

If you want to change the number of events to run over, go into `codes/functions.py`, scroll down to `run_ddsim` function, and edit `--numberOfEvents`

### Steps:
1. Run
    ```bash
    singularity run --nv --bind /cvmfs,/local /local/d1/badea/mu+mu-/mucoll-deploy.sif
    source /opt/setup_mucoll.sh
    ```
2. Make sure your samples are in `/local/d1/mu+mu-/samples` in the form mass_lifetime.hepmc (mass in GeV, lifetime in ns)
3. Navigate to the `sim_Hbb` directory
4. Run the following command:

   ```bash
   python multi_sim.py input_files -o output_directory -j number_of_cores
   ```
   Replace input_files with as many input files as you want in the form mass_lifetime.hepmc

   If -o is not specified, the default output directory is ../sim_Hbb
   
   If -j is not specified, the default is 2 cores
   
   Example command:
   ```bash
   python multi_sim.py 134_1.hepmc 300_10.hepmc -o /local/d1/lrozanov/mucoll-tutorial-2023/sim_Hbb -j 2
   ```

## Instructions for running multi_digi_reco.py:

### Steps:
1. Run
    ```bash
    singularity run --nv --bind /cvmfs,/local /local/d1/badea/mu+mu-/mucoll-deploy.sif
    source /opt/setup_mucoll.sh
    ```
2. Edit lines 53-58 to prepend the filepaths for whatever directories your files are in
3. Run the following command:

   ```bash
   python multi_digi_reco.py input_files -r reco -b bib -o output_directory -j number_of_cores
   ```
   Replace input_files with as many input files as you want in the form mass_lifetime.hepmc

   If -r is not included, the default is false

   If -b is not included, the default is false

   If -o is not specified, the default output directory is None
   
   If -j is not specified, the default is 1 cores
   
   Example command for digi w/ bib overlaid:
   ```bash
   python multi_digi_reco.py 134_1_sim.slcio 300_10_sim.slcio -b -o /local/d1/lrozanov/mucoll-tutorial-2023/digi_Hbb_bib -j 2
   ```
   Example command for reco w/ bib overlaid:
   ```bash
   python multi_digi_reco.py 134_1_digi_bib.slcio 300_10_digi_bib.slcio -r -b -o /local/d1/lrozanov/mucoll-tutorial-2023/reco_Hbb_bib -j 2
   ```