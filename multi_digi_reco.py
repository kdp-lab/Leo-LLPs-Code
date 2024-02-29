import argparse
import multiprocessing
import os
import subprocess
import re

# Define the run_marlin function within the script
def run_marlin(input_file, bib, reco, output_directory):
    base_dir = "/local/d1/lrozanov/mucoll-tutorial-2023"
    if reco:
        steering_file = f"{base_dir}/mucoll-benchmarks/reconstruction/marlin/reco_steer.xml"
        default_output_dir = f"{base_dir}/reco_Hbb{'_bib' if bib else ''}"
    else:
        steering_file = f"{base_dir}/mucoll-benchmarks/digitisation/marlin/digi_steer.xml"
        default_output_dir = f"{base_dir}/digi_Hbb{'_bib' if bib else ''}"

    if output_directory is None:
        output_directory = default_output_dir
    
    pattern = r'(_sim|_digi)\.slcio$'
    # Extract the base filename without the path
    base_filename = os.path.basename(input_file)
    # Use re.sub to remove '_sim.slcio' or '_digi.slcio' from the end of the filename
    input_filename = re.sub(pattern, '', base_filename)    

    overlay_option = "Test" if bib else "None"
    bib_suffix = "_bib" if bib else ""
    task_suffix = "reco" if reco else "digi"
    output_file_all = os.path.join(output_directory, f"{input_filename}_{task_suffix}{bib_suffix}.slcio")
    output_file_light = os.path.join(output_directory, f"{input_filename}_{task_suffix}{bib_suffix}_light.slcio")

    command = [
        "Marlin",
        steering_file,
        "--global.LCIOInputFiles="+input_file,
        "--DD4hep.DD4hepXMLFile="+os.getenv("MUCOLL_GEO"),
        "--Config.Overlay="+overlay_option,
        "--LCIOWriter_all.LCIOOutputFile="+output_file_all,
        "--LCIOWriter_light.LCIOOutputFile="+output_file_light
    ]
    subprocess.run(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Marlin tasks in parallel.")
    parser.add_argument("input_files", nargs="+", help="List of input files for processing.")
    parser.add_argument("-b", "--bib", action='store_true', help="Use the BIB overlay.", default=False)
    parser.add_argument("-r", "--reco", action='store_true', help="Run reconstruction instead of digitisation.", default=False)
    parser.add_argument("-o", "--output_directory", help="Output directory for task results.", default=None)
    parser.add_argument("-j", "--ncpu", help="Number of CPU cores to use.", type=int, default=1)
    args = parser.parse_args()

    # Prepend "/local/d1/lrozanov/mucoll-tutorial-2023/sim_Hbb/' to each input file path
    input_files = [f"/local/d1/lrozanov/mucoll-tutorial-2023/sim_Hbb/{input_file}" for input_file in args.input_files]

    # Use multiprocessing to parallelize task execution
    with multiprocessing.Pool(args.ncpu) as pool:
        pool.starmap(run_marlin, [(input_file, args.bib, args.reco, args.output_directory) for input_file in input_files])
