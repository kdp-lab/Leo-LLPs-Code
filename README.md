# Leo-LLPs-Code

This repository contains code for running multi-processing over sim, digi, and reco for studying long-lived particles (LLPs). 

It was necessary for me to rebuild `lcgeo` and `DD4HEP` after a new version was released with a fix for [this](https://github.com/AIDASoft/DD4hep/pull/1260) issue - TLDR: Geant4 could not read Hepmc status codes; now the status codes can be specified in the steering file.

When building `lcgeo` it was necessary to replace all `dd4hep::long64` with `long long int`.


### Prerequisites:
- Clone both [`MuC-tutorial`](https://github.com/MuonColliderSoft/MuC-Tutorial) into this directory.

Run the following commands for sim:

```bash
source /local/d1/brosser/setup.sh
source /local/d1/brosser/DD4hep/bin/thisdd4hep.sh
source /local/d1/brosser/lcgeo/bin/thislcgeo.sh
```

And the following commands for digi/reco:

```bash
singularity run --nv --bind /cvmfs,/local /local/d1/badea/mu+mu-/mucoll-deploy.sif
source /opt/setup_mucoll.sh
```

### Tips:
- If you're running over sim or digi_bib/reco_bib over many events, use `nohup` or `screen` to create a process that is separate from your connection to the server so that it can keep running in the background. 

E.g. `nohup my_python_script.py > output.log 2>&1 &`

##
 Instructions for running multi_sim.py:

### Steps:
1. Make sure your samples are in `/local/d1/mu+mu-/samples` in the form mass_lifetime.hepmc (mass in GeV, lifetime in ns)
2. Run the following command:

   ```bash
   python multi_sim.py input_files -o output_directory -j number_of_cores
   ```
   Replace input_files with as many input files as you want in the form mass_lifetime.hepmc

   You can use `-n` or `--number_of_events` to specify the number of events to simulate. The default is all events

   If `-o` is not specified, the default output directory is `/local/d1/mu+mu-/sim`
   
   If `-j` is not specified, the default is 2 cores
   
   Example command:
   ```bash
   python multi_sim.py 134_1.hepmc 300_10.hepmc -o /local/d1/lrozanov/mucoll-tutorial-2023/sim_Hbb -j 2
   ```

## Instructions for running multi_digi_reco.py:

### Steps:
(NB: You can be in any directory)
1. Make sure your input files are in `/local/d1/mu+mu-/...` `sim`, `digi`, `digi_bib` for whichever process you're running
2. Run the following command:

   ```bash
   python multi_digi_reco.py input_files -r reco -b bib -o output_directory -j number_of_cores
   ```
   Replace input_files with as many input files as you want in the form mass_lifetime.hepmc

   If `-r` is not included, the default is false

   If `-b` is not included, the default is false

   If `-o` is not specified, the default output directory is `/local/d1/mu+mu/digi` for digi, `/local/d1/mu+mu-/digi_bib` for digi with bib, `/local/d1/mu+mu-/reco` for reco, and `/local/d1/mu+mu-/reco_bib` for reco with bib.
   
   If `-j` is not specified, the default is 1 core
   
   Example command for digi w/ bib overlaid:
   ```bash
   python multi_digi_reco.py 134_1_sim.slcio 300_10_sim.slcio -b -o /local/d1/lrozanov/mucoll-tutorial-2023/digi_Hbb_bib -j 2
   ```
   Example command for reco w/ bib overlaid:
   ```bash
   python multi_digi_reco.py 134_1_digi_bib.slcio 300_10_digi_bib.slcio -r -b -o /local/d1/lrozanov/mucoll-tutorial-2023/reco_Hbb_bib -j 2
   ```
