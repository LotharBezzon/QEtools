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