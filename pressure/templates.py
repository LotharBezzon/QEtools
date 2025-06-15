import os

QE_VERSION = '7.4.1'
PSEUDO_DIR = os.path.join(f'qe-{QE_VERSION}', 'pseudo')
OUTDIR_BASE = "./tmp_output_python/"

NATOMS = 2
NTYPES = 1

CELL_PARAMS = "  0.0000 2.8369 2.8369\n  2.8369 0.0000 2.8369\n  2.8369 2.8369 0.0000\n"
ATOMIC_POSITIONS = ["1.4184 4.2553 1.4184", "0.0000 0.0000 2.8369"]
DEGAUSS = 0.02

# Converged parameters (replace with your actual converged values)
ECUTWFC = 50.0
ECUTRHO = 200.0
K_POINTS_GRID = "6 6 6 1 1 1"
PSEUDOPOTENTIAL_FILE = "Ge.pbe-kjpaw.UPF"

Ge_I_params = {
    'name': 'Ge-I',
    'pseudo_dir': PSEUDO_DIR,
    'outdir_base': OUTDIR_BASE,
    'natoms': NATOMS,
    'ntypes': NTYPES,
    'ecutwfc': ECUTWFC,
    'ecutrho': ECUTRHO,
    'cell_params': CELL_PARAMS,
    'atomic_positions': ATOMIC_POSITIONS,
    'k_points_grid': K_POINTS_GRID,
    'degauss': 0.0,
    'metal': False
}

Ge_II_params = {
    'name': 'Ge-II',
    'pseudo_dir': PSEUDO_DIR,
    'outdir_base': OUTDIR_BASE,
    'natoms': 2,
    'ntypes': 1,
    'ecutwfc': ECUTWFC,
    'ecutrho': ECUTRHO,
    'cell_params': "5.0000 0.0000 0.0000\n0.0000 5.0000 0.0000\n0.0000 0.0000 2.7500\n",
    'atomic_positions': ["0.00 0.00 0.25", "0.50 0.50 0.75"],
    'k_points_grid': K_POINTS_GRID,
    'degauss': 0.02,
    'metal': True
}


def prepare_template(name, pseudo_dir, outdir_base, natoms, ntypes, ecutwfc, ecutrho, cell_params, atomic_positions, k_points_grid, metal, degauss=0.0, smearing='gaussian'):
    """
    Prepare the QE input template with the given parameters.
    """
    template = f"""
&CONTROL
  calculation = 'vc-relax'
  prefix = '{name}'
  outdir = '{outdir_base}'
  pseudo_dir = '{pseudo_dir}'
  tprnfor = .true.
  tstress = .true.
  etot_conv_thr = 1.0d-4
  forc_conv_thr = 1.0d-3
/
&SYSTEM
  ibrav = 0
  nat = {natoms}
  ntyp = {ntypes}
  ecutwfc = {ecutwfc}
  ecutrho = {ecutrho}
"""
    if metal:
        template += f"""
  occupation = 'smearing'
  smearing = {smearing}
  degauss = {degauss}
"""
    template += f"""
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
CELL_PARAMETERS {{angstrom}}
{cell_params}
ATOMIC_SPECIES
  Ge  72.64  {PSEUDOPOTENTIAL_FILE}
ATOMIC_POSITIONS {{angstrom}}
  Ge  {atomic_positions[0]}
  Ge  {atomic_positions[1]}
K_POINTS {{automatic}}
  {k_points_grid}
"""
    return template

# --- QE Input Template ---
# Use f-strings for easy variable substitution later
Ge_I = f"""
&CONTROL
  calculation = 'vc-relax'
  prefix = 'Ge-I'
  outdir = '{OUTDIR_BASE}'
  pseudo_dir = '{PSEUDO_DIR}'
  tprnfor = .true.
  tstress = .true.
  etot_conv_thr = 1.0d-4
  forc_conv_thr = 1.0d-3
/
&SYSTEM
  ibrav = 0
  nat = {NATOMS}
  ntyp = {NTYPES}
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
CELL_PARAMETERS {{angstrom}}
{CELL_PARAMS}
ATOMIC_SPECIES
  Ge  72.64  {PSEUDOPOTENTIAL_FILE}
ATOMIC_POSITIONS {{angstrom}}
  Ge  {ATOMIC_POSITIONS[0]}
  Ge  {ATOMIC_POSITIONS[1]}
K_POINTS {{automatic}}
  {K_POINTS_GRID}
"""

Ge_II = f"""
&CONTROL
  calculation = 'vc-relax'
  prefix = 'Ge-II'
  outdir = '{OUTDIR_BASE}'
  pseudo_dir = '{PSEUDO_DIR}'
  tprnfor = .true.
  tstress = .true.
  etot_conv_thr = 1.0d-4
  forc_conv_thr = 1.0d-3
/
&SYSTEM
  ibrav = 0
  nat = 2
  ntyp = 1
  ecutwfc = {ECUTWFC}
  ecutrho = {ECUTRHO}
  occupation = 'smearing'
  smearing = 'gaussian'
  degauss = 0.02
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
  Ge  72.64  {PSEUDOPOTENTIAL_FILE}
CELL_PARAMETERS (angstrom)
5.000000000   0.000000000   0.000000000
0.000000000   5.000000000   0.000000000
0.000000000   0.000000000   2.750000000
! Initial guess for Ge-II lattice parameters (a ~ 5.0 Å, c ~ 2.75 Å for c/a ~ 0.55)
! These are highly dependent on the target pressure.
ATOMIC_POSITIONS (crystal)
Ge  0.00  0.00  0.25
Ge  0.50  0.50  0.75
K_POINTS (automatic)
  {K_POINTS_GRID}
"""