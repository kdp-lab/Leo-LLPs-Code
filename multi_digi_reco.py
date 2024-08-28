import argparse
import multiprocessing
import os
import subprocess
import re

# Define the run_marlin function within the script
def run_marlin(input_file, bib, reco, output_directory, number_of_events):
    base_dir = "/local/d1/lrozanov/mucoll-tutorial-2023"
    if reco:
        steering_file = f"{base_dir}/mucoll-benchmarks-LLPs/reconstruction/marlin/reco_steer.xml" # HAVE TO EXPORT FIRST: cd /local/d1/mu+mu-/ and then export k4run_DLL=$(realpath libMyBIBUtils.so):${MARLIN_DLL}
        default_output_dir = f"{base_dir}/reco_Hbb{'_bib' if bib else ''}"
    else:
        steering_file = f"{base_dir}/mucoll-benchmarks-LLPs/digitisation/marlin/digi_steer.xml" # https://github.com/mlarson02/mucoll-benchmarks-LLPs/tree/main
        default_output_dir = f"{base_dir}/digi_Hbb{'_bib' if bib else ''}"

    if output_directory is None:
        output_directory = default_output_dir
    
    pattern = r'(_sim|_digi|_digi_bib)\.slcio$'
    # Extract the base filename without the path
    base_filename = os.path.basename(input_file)
    # Use re.sub to remove '_sim.slcio' or '_digi.slcio' from the end of the filename
    input_filename = re.sub(pattern, '', base_filename)    

    bib_suffix = "_bib" if bib else ""
    task_suffix = "reco" if reco else "digi"
    output_file_all = os.path.join(output_directory, f"{input_filename}_{task_suffix}{bib_suffix}.slcio")
    output_file_light = os.path.join(output_directory, f"{input_filename}_{task_suffix}{bib_suffix}_light.slcio")
    output_file_root = os.path.join(output_directory, f"{input_filename}_{task_suffix}{bib_suffix}")

    command = [
        "Marlin",
        steering_file,
        "--global.LCIOInputFiles="+input_file,
        "--DD4hep.DD4hepXMLFile="+os.getenv("MUCOLL_GEO"),
        "--LCIOWriter_all.LCIOOutputFile="+output_file_all,
        "--LCIOWriter_light.LCIOOutputFile="+output_file_light,
        "--AIDA.FileName="+output_file_root
    ]
    if bib:
        command += ["--Config.Overlay=Test"]
    if number_of_events > 0:
        command += ["--global.MaxRecordNumber="+f"{number_of_events+1}"] 
        
    subprocess.run(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Marlin tasks in parallel.")
    parser.add_argument("input_files", nargs="+", help="List of input files for processing.")
    parser.add_argument("-n", "--number_of_events", help="Number of events to simulate.", type=int, default=-1)
    parser.add_argument("-b", "--bib", action='store_true', help="Use the BIB overlay.", default=False)
    parser.add_argument("-r", "--reco", action='store_true', help="Run reconstruction instead of digitisation.", default=False)
    parser.add_argument("-o", "--output_directory", help="Output directory for task results.", default= "/local/d1/mu+mu-/digi_v2")
    parser.add_argument("-i", "--input_directory", help="Input directory for sim/digi files if not in /local/d1/mu+mu-/sim or /digi_v2.", default= "/local/d1/mu+mu-/sim")
    parser.add_argument("-j", "--ncpu", help="Number of CPU cores to use.", type=int, default=1)
    args = parser.parse_args()


    if args.reco and args.bib:
        base_directory = "/local/d1/mu+mu-/digi_bib_v2"
        output_directory = "/local/d1/mu+mu-/reco_bib_v2"
    elif args.reco:
        base_directory = "/local/d1/mu+mu-/digi_v2"
        output_directory = "/local/d1/mu+mu-/reco_v2"
    else:
        base_directory = "/local/d1/mu+mu-/sim"
        output_directory = "/local/d1/mu+mu-/digi_v2"

    if args.output_directory != "/local/d1/mu+mu-/digi_v2":
        output_directory = args.output_directory
        
    if args.input_directory != "/local/d1/mu+mu-/sim" and args.input_directory != "/local/d1/mu+mu-/digi_v2":
        base_directory = args.input_directory

    # Prepend the above to each input file path
    input_files = [f"{base_directory}/{input_file}" for input_file in args.input_files]

    # Use multiprocessing to parallelize task execution
    with multiprocessing.Pool(args.ncpu) as pool:
        pool.starmap(run_marlin, [(input_file, args.bib, args.reco, output_directory, args.number_of_events) for input_file in input_files])
