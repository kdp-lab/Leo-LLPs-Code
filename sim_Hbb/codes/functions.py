import os
import subprocess


def extract_params(filename):
    # Extract mass and lifetime from the filename
    mass_lifetime = os.path.splitext(os.path.basename(filename))[0]  # Remove file extension
    mass_str, lifetime_str = mass_lifetime.split('_')
    
    # Convert to float
    mass = float(mass_str)  # Assuming mass is in GeV
    lifetime = float(lifetime_str)  # Assuming lifetime is in ns
    
    # Calculate width in GeV
    hbar = 6.582219514e-16  # Planck constant in eV·s
    width_ev = hbar / (lifetime* 1e-9)  # Width in eV, lifetime in s
    width = width_ev * (10 ** (-9))  # Width in GeV

    # Convert lifetime to mm
    lifetime = lifetime * 299.792458  # Conversion from ns to mm

    return mass, width, lifetime


def generate_modified_tbl(mass, lifetime, width):
    input_tbl_file = "/local/d1/lrozanov/mucoll-tutorial-2023/sim_Hbb/tbl_files/example.tbl"
    output_tbl_file = f"tbl_files/{mass:.0f}GeV_{lifetime:.0f}mm.tbl"
    
    # Check if the output directory exists
    output_dir = os.path.dirname(output_tbl_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create directory if it doesn't exist

    # Check if the output file exists
    if os.path.exists(output_tbl_file):
        print(f"Output file '{output_tbl_file}' already exists. Skipping generation.")
        return output_tbl_file
    
    with open(input_tbl_file, 'r') as infile:
        lines = infile.readlines()

    # Find the index of the line containing '1000015' or '2000015'
    index_1000015 = next((i for i, line in enumerate(lines) if '1000015' in line), None)
    index_2000015 = next((i for i, line in enumerate(lines) if '2000015' in line), None)

    # Modify the line containing '1000015' or '2000015' and the subsequent line
    if index_1000015 is not None:
        if index_1000015 < len(lines) - 1:
            parts = lines[index_1000015].split()
            parts[3] = f"{mass:.5f}"  # Mass
            parts[4] = f"{width:.5e}"  # Width
            parts[5] = f"{lifetime:.5f}"  # Lifetime (converted to mm)
            updated_line = "{:>11} {:<22} {:<5} {:<13} {:<12} {:<12}\n".format(*parts)
            lines[index_1000015] = updated_line

            # Modify the line below
            parts_next = lines[index_1000015 + 1].split()
            parts_next[3] = f"{mass:.5f}"  # Mass
            parts_next[4] = f"{width:.5e}"  # Width
            parts_next[5] = f"{lifetime:.5f}"  # Lifetime (converted to mm)
            updated_line = "{:>11} {:<23} {:<4} {:<13} {:<12} {:<12}\n".format(*parts_next)
            lines[index_1000015 + 1] = updated_line

    if index_2000015 is not None:
        if index_2000015 < len(lines) - 1:
            parts = lines[index_2000015].split()
            parts[3] = f"{mass:.5f}"  # Mass
            parts[4] = f"{width:.5e}"  # Width
            parts[5] = f"{lifetime:.5f}"  # Lifetime (converted to mm)
            updated_line = "{:>11} {:<22} {:<5} {:<13} {:<12} {:<12}\n".format(*parts)
            lines[index_2000015] = updated_line

            # Modify the line below
            parts_next = lines[index_2000015 + 1].split()
            parts_next[3] = f"{mass:.5f}"  # Mass
            parts_next[4] = f"{width:.5e}"  # Width
            parts_next[5] = f"{lifetime:.5f}"  # Lifetime (converted to mm)
            updated_line = "{:>11} {:<23} {:<4} {:<13} {:<12} {:<12}\n".format(*parts_next)
            lines[index_2000015 + 1] = updated_line

    # Write modified content to output file
    with open(output_tbl_file, 'w') as outfile:
        outfile.writelines(lines)

    return output_tbl_file

def run_ddsim(input_file, output_directory, tbl_file, number_of_events):
    # Extracting the filename without extension
    input_filename = os.path.splitext(os.path.basename(input_file))[0]
    if f"{input_filename}_sim.slcio" in os.listdir(output_directory):
        print(f"Overriding existing file: {input_filename}_sim.slcio")
        subprocess.run(["rm", f"{output_directory}/{input_filename}_sim.slcio", "-f"])
    # DDSIM command
    command = [
        "ddsim",
        "--steeringFile",
        "/local/d1/lrozanov/mucoll-tutorial-2023/mucoll-benchmarks/simulation/ilcsoft/steer_baseline.py",
        "--inputFile",
        input_file,
        "--physics.pdgfile",
        f"{tbl_file}",
        "--outputFile",
        f"{output_directory}/{input_filename}_sim.slcio"
        # ,"--dumpParameter", "--dump" # Uncomment to print out the parameters in steering file
    ]
    if number_of_events > 0:
        command += ["--numberOfEvents", f"{number_of_events}"]
        
    subprocess.run(command)