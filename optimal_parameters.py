import os
import subprocess
import re
import matplotlib.pyplot as plt
import numpy as np
from templates import *
import argparse

parser = argparse.ArgumentParser(description="Run QE scf calculations to find optimal ecutwfc and k-points parameters.")
parser.add_argument('-t', '--template', type=str)
args = parser.parse_args()

if args.template+'_params' in globals():
    QE_INPUT_TEMPLATE = globals()[args.template+'_params']
    del QE_INPUT_TEMPLATE['ecutwfc', 'ecutrho', 'k_points_grid']    # These will be tested later
else:
    raise ValueError(f"Template '{args.template}' not found. Add it to templates.py.")

def generate_input_file(input_template, filename, ecutwfc, k_points):
    """
    Generates a Quantum ESPRESSO input file for an SCF calculation.
    """
    calculation = 'scf'
    input_content = prepare_input(**input_template, calculation=calculation, ecutwfc=ecutwfc, ecutrho=4*ecutwfc, k_points_grid=k_points)
    with open(filename, 'w') as f:
        f.write(input_content)

def run_pwscf(input_file, output_file, pw_x_path='pw.x', n_cores=2):
    """
    Runs the pw.x executable.

    Args:
        input_file (str): Path to the Quantum ESPRESSO input file.
        output_file (str): Path where the output will be written.
        pw_x_path (str): Path to the pw.x executable.

    Returns:
        bool: True if the command executed successfully, False otherwise.
    """
    print(f"Running {pw_x_path} with {input_file}...")
    try:
        command = f"mpirun -np {n_cores} {pw_x_path} -inp {input_file}"
        with open(output_file, 'w') as outfile:
            subprocess.run(command.split(), stdout=outfile, stderr=subprocess.PIPE, check=True)
        print(f"Successfully ran {input_file}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running pw.x: {e}")
        print(f"Stderr: {e.stderr.decode()}")
        return False
    except FileNotFoundError:
        print(f"Error: pw.x executable not found at '{pw_x_path}'. "
              "Please ensure it's in your PATH or provide the full path.")
        return False

def parse_total_energy(output_file):
    """
    Parses the Quantum ESPRESSO output file to extract the final total energy.

    Args:
        output_file (str): Path to the Quantum ESPRESSO output file.

    Returns:
        float or None: Total energy in Rydberg, or None if not found.
    """
    with open(output_file, 'r') as f:
        content = f.read()
        # Regex to find the last occurrence of '!    total energy'
        match = re.search(r'!\s+total energy\s+=\s+([-+]?\d+\.\d+)\s+Ry', content)
        if match:
            return float(match.group(1))
        return None

def plot_convergence(x_values, energies, xlabel, title, output_filename):
    """
    Plots the convergence curve.

    Args:
        x_values (list): List of x-axis values (e.g., cutoff or k-points).
        energies (list): List of corresponding total energies.
        xlabel (str): Label for the x-axis.
        title (str): Title of the plot.
        output_filename (str): Filename to save the plot.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, energies, marker='o', linestyle='-', color='blue')
    plt.xlabel(xlabel)
    plt.ylabel('Total Energy (Ry)')
    plt.title(title)
    plt.grid(True)
    plt.axhline(y=min(energies), color='red', linestyle='--', label='Min Energy')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_filename)
    print(f"Convergence plot saved to {output_filename}")
    # plt.show() # Uncomment to display plot immediately

def find_converged_value(x_values, energies, tolerance=1e-4):
    """
    Finds the converged x-value based on energy change.

    Args:
        x_values (list): List of x-axis values.
        energies (list): List of corresponding total energies.
        tolerance (float): Maximum allowed energy difference between consecutive points (in Ry).

    Returns:
        tuple: (converged_x_value, converged_energy) or (None, None) if not converged.
    """
    if len(energies) < 2:
        return None, None

    for i in range(len(energies) - 1, 0, -1):
        # Check from the end of the list
        if abs(energies[i] - energies[i-1]) < tolerance:
            # Check if the next point also meets the criteria (to avoid false positives on small dips)
            if i + 1 < len(energies) and abs(energies[i+1] - energies[i]) < tolerance:
                return x_values[i], energies[i]
            elif i + 1 == len(energies): # If it's the last point, and it converged
                 return x_values[i], energies[i]
    return None, None


def main():
    # --- Configuration ---
    PW_X_PATH = './qe-7.4.1/bin/pw.x'  # Make sure 'pw.x' is in your system's PATH or provide full path
    PSEUDO_DIR = './qe-7.4.1/pseudo/' # Directory where your .UPF pseudopotential files are located
    SYSTEM_NAME = args.template

    # Create a temporary directory for output files
    os.makedirs('./tmp', exist_ok=True)
    os.makedirs(PSEUDO_DIR, exist_ok=True) # Ensure pseudo_dir exists, though actual UPF file needs to be there

    print("Starting Quantum ESPRESSO Convergence Tests...")

    # --- 1. Ecutwfc Convergence Test ---
    print("\n--- Ecutwfc Convergence ---")
    # Define ecutwfc values to test (adjust based on your pseudopotential type)
    # For norm-conserving, start lower (e.g., 20-40 Ry). For ultrasoft/PAW, can be higher.
    ecutwfc_values = list(range(20, 81, 5)) # Example: 20, 25, ..., 80 Ry
    fixed_k_points = "4 4 4 1 1 1" # A reasonable starting k-point grid for many bulk systems

    ecutwfc_energies = []
    ecutwfc_tested = []

    for ecut in ecutwfc_values:
        input_file = f"tmp/{SYSTEM_NAME}_ecut_{ecut}.in"
        output_file = f"tmp/{SYSTEM_NAME}_ecut_{ecut}.out"

        generate_input_file(QE_INPUT_TEMPLATE, input_file, ecut, fixed_k_points)

        if run_pwscf(input_file, output_file, PW_X_PATH):
            energy = parse_total_energy(output_file)
            if energy is not None:
                ecutwfc_energies.append(energy)
                ecutwfc_tested.append(ecut)
                print(f"Ecutwfc: {ecut} Ry, Energy: {energy:.6f} Ry")
            else:
                print(f"Warning: Could not parse energy for ecutwfc = {ecut}")
        else:
            print(f"Skipping ecutwfc = {ecut} due to pw.x error.")

    if not ecutwfc_tested:
        print("Ecutwfc convergence test failed: no data collected.")
        return

    plot_convergence(ecutwfc_tested, ecutwfc_energies, 'Ecutwfc (Ry)',
                     f'{SYSTEM_NAME} Ecutwfc Convergence', f'{SYSTEM_NAME}_ecutwfc_convergence.png')

    converged_ecutwfc, _ = find_converged_value(ecutwfc_tested, ecutwfc_energies, tolerance=5e-4) # Ry tolerance for convergence
    if converged_ecutwfc:
        print(f"\nRecommended Ecutwfc: {converged_ecutwfc} Ry (converged to {5e-4} Ry)")
    else:
        print(f"\nCould not find a clear converged Ecutwfc. Review the plot and adjust range/tolerance if needed.")
        converged_ecutwfc = ecutwfc_tested[-1] # Fallback to highest tested if no convergence detected

    # --- 2. K-points Convergence Test ---
    print("\n--- K-points Convergence ---")
    # Define k-point mesh values (e.g., for nxnxn 1 1 1 grid)
    # Start from a coarse grid and increase density
    k_mesh_values = list(range(2, 11, 1)) # Example: 2x2x2 to 10x10x10 grids

    k_points_energies = []
    k_points_tested = []

    if converged_ecutwfc is None:
        print("Cannot proceed with k-points convergence without a converged ecutwfc.")
        return

    for n_k in k_mesh_values:
        input_file = f"tmp/{SYSTEM_NAME}_kpt_{n_k}.in"
        output_file = f"tmp/{SYSTEM_NAME}_kpt_{n_k}.out"
        current_k_points = f"{n_k} {n_k} {n_k} 1 1 1" # Using shifted Monkhorst-Pack grid

        generate_input_file(QE_INPUT_TEMPLATE, input_file, converged_ecutwfc, current_k_points)

        if run_pwscf(input_file, output_file, PW_X_PATH):
            energy = parse_total_energy(output_file)
            if energy is not None:
                k_points_energies.append(energy)
                k_points_tested.append(n_k)
                print(f"K-mesh: {n_k}x{n_k}x{n_k}, Energy: {energy:.6f} Ry")
            else:
                print(f"Warning: Could not parse energy for k-mesh = {n_k}")
        else:
            print(f"Skipping k-mesh = {n_k} due to pw.x error.")

    if not k_points_tested:
        print("K-points convergence test failed: no data collected.")
        return

    plot_convergence(k_points_tested, k_points_energies, 'K-point Mesh Size (n)',
                     f'{SYSTEM_NAME} K-points Convergence (using Ecutwfc={converged_ecutwfc} Ry)',
                     f'{SYSTEM_NAME}_kpoints_convergence.png')

    converged_k_mesh_n, _ = find_converged_value(k_points_tested, k_points_energies, tolerance=1e-4) # Ry tolerance for convergence
    if converged_k_mesh_n:
        print(f"\nRecommended K-point Mesh (n): {converged_k_mesh_n}x{converged_k_mesh_n}x{converged_k_mesh_n} (converged to {1e-4} Ry)")
        print(f"Full K_POINTS line: {converged_k_mesh_n} {converged_k_mesh_n} {converged_k_mesh_n} 1 1 1")
    else:
        print(f"\nCould not find a clear converged K-point mesh. Review the plot and adjust range/tolerance if needed.")
        converged_k_mesh_n = k_points_tested[-1] # Fallback to highest tested if no convergence detected
        print(f"Using highest tested K-point mesh (n): {converged_k_mesh_n}x{converged_k_mesh_n}x{converged_k_mesh_n}")


    print("\n--- Summary of Recommended Parameters ---")
    print(f"Converged Ecutwfc: {converged_ecutwfc} Ry")
    print(f"Converged K-point Mesh: {converged_k_mesh_n}x{converged_k_mesh_n}x{converged_k_mesh_n} 1 1 1")
    print("\nRemember to manually verify the convergence plots for visual confirmation.")

if __name__ == '__main__':
    main()
