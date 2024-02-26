# Leo-LLPs-Code

This repository contains code for running multi-simulations using DDSIM for studying long-lived particles (LLPs).

## Instructions for running multi_sim.py:

### Prerequisites:
- Clone both [`MuC-tutorial`](https://github.com/MuonColliderSoft/MuC-Tutorial) and [`mucoll-benchmarks`](https://github.com/MuonColliderSoft/mucoll-benchmarks/tree/main) into this directory


### Steps:
1. Run
    ```bash
    singularity run --nv --bind /cvmfs,/local /local/d1/badea/mu+mu-/mucoll-deploy.sif
    source /opt/setup_mucoll.sh
    ```
2. Make sure your samples are in `/local/d1/mu+mu-/samples`
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

