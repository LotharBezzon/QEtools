import os

QE_VERSION = '7.4.1'
PSEUDO_DIR = os.path.join(f'qe-{QE_VERSION}', 'pseudo')
OUTDIR_BASE = "./tmp_output_python/"

# Initial cell parameter for Ge in Bohr (approximate, vc-relax will optimize)
# This is celldm(1) in QE. For Ge, it's around 10.7 Bohr at 0 GPa.
INITIAL_CELLDM1 = 10.7 

# Converged parameters (replace with your actual converged values)
ECUTWFC = 40.0
ECUTRHO = 160.0
K_POINTS_GRID = "6 6 6 0 0 0"
PSEUDOPOTENTIAL_FILE = "Ge.pbe-kjpaw.UPF"

# --- QE Input Template ---
# Use f-strings for easy variable substitution later
Ge_cd = f"""
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
  press_conv_thr = 0.5  ! kbar
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