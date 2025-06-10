import os
import subprocess
import re # For regular expressions to parse output

# --- Configuration ---
QE_BIN = "/path/to/your/qe/bin/pw.x" 
PSEUDO_DIR = "./pseudos/"
OUTDIR_BASE = "./tmp_output_python/"
NUM_CORES = 4 # Adjust based on your laptop's CPU

# Pressures in kbar (1 GPa = 10 kbar)
PRESSURES_KBAR = [-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100]

# Initial cell parameter for Ge in Bohr (approximate, vc-relax will optimize)
# This is celldm(1) in QE. For Ge, it's around 10.7 Bohr at 0 GPa.
INITIAL_CELLDM1 = 10.7 

# Converged parameters (replace with your actual converged values)
ECUTWFC = 40.0
ECUTRHO = 160.0
K_POINTS_GRID = "6 6 6 0 0 0"
PSEUDOPOTENTIAL_FILE = "Ge.pbe-spn-kjpaw_psl.1.0.0.UPF"

# --- QE Input Template ---
# Use f-strings for easy variable substitution later
QE_INPUT_TEMPLATE = f"""
&CONTROL
  calculation = 'vc-relax'
  prefix = 'Ge'
  outdir = '{OUTDIR_BASE}'
  pseudo_dir = '{PSEUDO_DIR}'
  tprnfor = .true.
  tstress = .true.
  etot_conv_thr = 1.0d-4
  forc_conv_thr = 1.0d-3
/
&SYSTEM
  ibrav = 2,
  celldm(1) = {INITIAL_CELLDM1}
  nat = 2
  ntyp = 1
  ecutwfc = {ECUTWFC}
  ecutrho = {ECUTRHO}
/
&ELECTRONS
  conv_thr = 1.0d-8
  mixing_beta = 0.7
/
&IONS
  ion_dynamics = 'bfgs'
/
&CELL
  cell_dynamics = 'bfgs'
  press = {{pressure_val}} ! Placeholder for pressure
  press_conv_thr = 0.5
  cell_dofree = 'all'
/
ATOMIC_SPECIES
  Ge  72.63  {PSEUDOPOTENTIAL_FILE}
ATOMIC_POSITIONS (crystal)
  Ge  0.0  0.0  0.0
  Ge  0.25 0.25 0.25
K_POINTS (automatic)
  {K_POINTS_GRID}
"""

# --- Main Script Logic ---
def run_qe_relaxation(pressure_kbar):
    """
    Generates QE input, runs vc-relax, and extracts key data.
    """
    run_dir = f"pressure_{pressure_kbar}kbar"
    os.makedirs(run_dir, exist_ok=True) # Create directory if it doesn't exist

    input_filename = os.path.join(run_dir, f"Ge_vc-relax_{pressure_kbar}kbar.in")
    output_filename = os.path.join(run_dir, f"Ge_vc-relax_{pressure_kbar}kbar.out")

    # Generate input file content with current pressure
    input_content = QE_INPUT_TEMPLATE.format(pressure_val=pressure_kbar)

    with open(input_filename, 'w') as f:
        f.write(input_content)

    print(f"Running relaxation for P = {pressure_kbar} kbar in {run_dir}...")

    # Run Quantum ESPRESSO using subprocess
    try:
        command = f"mpirun -np {NUM_CORES} {QE_BIN} -inp {input_filename}"
        with open(output_filename, 'w') as outfile:
            process = subprocess.run(command.split(), stdout=outfile, stderr=subprocess.PIPE, check=True)
        
        print(f"Calculation for P = {pressure_kbar} kbar finished.")
        
        # --- Post-processing: Extract data from output file ---
        final_enthalpy = None
        final_volume_au3 = None
        final_lattice_param_bohr = None
        
        with open(output_filename, 'r') as f:
            for line in f:
                if "Final enthalpy" in line:
                    match = re.search(r"Final enthalpy\s*=\s*([-\d.]+)\s*Ry", line)
                    if match:
                        final_enthalpy = float(match.group(1))
                elif "new unit-cell volume" in line:
                    match = re.search(r"new unit-cell volume\s*=\s*([-\d.]+)\s*a.u.\^3", line)
                    if match:
                        final_volume_au3 = float(match.group(1))
                elif "CELL_PARAMETERS" in line and "alat=" in line:
                    # This line is trickier if you're not explicitly tracking alat.
                    # For ibrav=2, the lattice constant is celldm(1) * 0.5.
                    # Or, easier, search for "Lattice constant =" line usually near the end
                    # For ibrav=2 (diamond cubic), the cubic lattice parameter 'a' is given by celldm(1)
                    # This value is usually printed at the end of a vc-relax run
                    pass # Will search for explicit lattice constant line
                elif "Lattice constant =" in line: # Often appears at the end for cubic systems
                    match = re.search(r"Lattice constant\s*=\s*([-\d.]+)\s*Bohr", line)
                    if match:
                        final_lattice_param_bohr = float(match.group(1))

        # If not found using "Lattice constant =", calculate from volume assuming cubic
        if final_lattice_param_bohr is None and final_volume_au3 is not None:
             # For cubic, V = a^3, so a = V^(1/3)
            final_lattice_param_bohr = final_volume_au3**(1/3)
            print(f"  (Calculated lattice param from volume: {final_lattice_param_bohr:.4f} Bohr)")

        return {
            'pressure_kbar': pressure_kbar,
            'enthalpy_Ry': final_enthalpy,
            'volume_au3': final_volume_au3,
            'lattice_param_bohr': final_lattice_param_bohr
        }

    except subprocess.CalledProcessError as e:
        print(f"Error running QE for P = {pressure_kbar} kbar: {e}")
        print(f"Stderr: {e.stderr.decode()}")
        return None
    except FileNotFoundError:
        print(f"Error: pw.x not found at {QE_BIN}. Please check the path.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None

if __name__ == "__main__":
    # Create main output directory if it doesn't exist
    os.makedirs(OUTDIR_BASE, exist_ok=True)
    os.makedirs(PSEUDO_DIR, exist_ok=True) # Make sure pseudo dir exists

    # Make sure your pseudopotential is in the PSEUDO_DIR
    if not os.path.exists(os.path.join(PSEUDO_DIR, PSEUDOPOTENTIAL_FILE)):
        print(f"Error: Pseudopotential file '{PSEUDOPOTENTIAL_FILE}' not found in '{PSEUDO_DIR}'.")
        print("Please download it and place it there.")
        exit(1)

    all_results = []
    for p in PRESSURES_KBAR:
        result = run_qe_relaxation(p)
        if result:
            all_results.append(result)

    # --- Print Summary (Optional) ---
    print("\n--- Summary of Results ---")
    print(f"{'Pressure (kbar)':<15} {'Enthalpy (Ry)':<18} {'Volume (a.u.^3)':<18} {'Lattice (Bohr)':<18}")
    print("-" * 70)
    for res in all_results:
        print(f"{res['pressure_kbar']:<15.1f} {res['enthalpy_Ry']:<18.8f} {res['volume_au3']:<18.8f} {res['lattice_param_bohr']:<18.8f}")

    # --- Further Analysis (e.g., plot Equation of State) ---
    # You would typically use libraries like NumPy and Matplotlib for this.
    # Convert results to arrays for plotting
    try:
        import numpy as np
        import matplotlib.pyplot as plt
        
        pressures = np.array([r['pressure_kbar'] for r in all_results])
        enthalpies = np.array([r['enthalpy_Ry'] for r in all_results])
        volumes = np.array([r['volume_au3'] for r in all_results])

        # Convert kbar to GPa for plotting if preferred (1 GPa = 10 kbar)
        pressures_gpa = pressures / 10.0

        plt.figure(figsize=(10, 6))
        plt.plot(pressures_gpa, enthalpies, 'o-', label='Enthalpy vs. Pressure')
        plt.xlabel('Pressure (GPa)')
        plt.ylabel('Enthalpy (Ry)')
        plt.title('Enthalpy-Pressure Curve for Germanium')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig('enthalpy_pressure_curve.png')
        plt.show()

        plt.figure(figsize=(10, 6))
        plt.plot(pressures_gpa, volumes, 'o-', label='Volume vs. Pressure')
        plt.xlabel('Pressure (GPa)')
        plt.ylabel('Volume (a.u.^3)')
        plt.title('Equation of State for Germanium')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig('volume_pressure_curve.png')
        plt.show()

    except ImportError:
        print("\nMatplotlib or NumPy not installed. Cannot generate plots.")
        print("Install with: pip install numpy matplotlib")

    print("\nScript finished.")