import subprocess
import multiprocessing
import argparse
import os
from codes import functions

def process_input_file(input_file, output_directory):
    # Extract mass, width, and lifetime from the input file
    mass, width, lifetime = functions.extract_params(input_file)
    # Generate modified .tbl file
    tbl_file = functions.generate_modified_tbl(mass, lifetime, width)
    # Run DDSIM simulation
    functions.run_ddsim(input_file, output_directory, tbl_file)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run multiple DDSIM simulations simultaneously.")
    parser.add_argument("input_files", nargs="+", help="List of input files for simulations.")
    parser.add_argument("-o", "--output_directory", help="Output directory for simulation results.", default="../sim_Hbb")
    parser.add_argument("-j", "--ncpu", help="Number of CPU cores to use.", type=int, default=1)
    args = parser.parse_args()
    
    # Prepend '../samples/' to each input file path
    input_files = [f"../samples/{input_file}" for input_file in args.input_files]

    # Create a pool of processes for parallel execution
    pool = multiprocessing.Pool(args.ncpu)

    # Process each input file
    pool.starmap(process_input_file, [(input_file, args.output_directory) for input_file in input_files])

    # Close the pool of processes
    pool.close()
    pool.join()

