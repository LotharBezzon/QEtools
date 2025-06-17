import os

QE_VERSION = '7.4.1'
PSEUDO_DIR = os.path.join(f'qe-{QE_VERSION}', 'pseudo')
OUTDIR_BASE = "./tmp_output_python/"

NATOMS = 2
NTYPES = 1

CELL_PARAMS = "  0.0000 2.8369 2.8369\n  2.8369 0.0000 2.8369\n  2.8369 2.8369 0.0000\n"
ATOMIC_POSITIONS = "  Ge  1.4184 4.2553 1.4184\n  Ge  0.0000 0.0000 2.8369"
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
    'metal': False,
    'atoms': "  Ge  72.64  Ge.pbe-kjpaw.UPF"
}

Ge_II_params = {
    'name': 'Ge-II',
    'pseudo_dir': PSEUDO_DIR,
    'outdir_base': OUTDIR_BASE,
    'natoms': 2,
    'ntypes': 1,
    'ecutwfc': ECUTWFC,
    'ecutrho': ECUTRHO,
    'cell_params': "  -2.5519 2.5519 1.4071\n  2.5519 -2.5519 1.4071\n  2.5519 2.5519 -1.4071\n",
    'atomic_positions': "  Ge  2.5519 2.5519 -0.0000\n  Ge  2.5519 -0.0000 0.7035",
    'k_points_grid': K_POINTS_GRID,
    'degauss': 0.02,
    'metal': True,
    'atoms': "  Ge  72.64  Ge.pbe-kjpaw.UPF"
}


def prepare_input(name, calculation, pseudo_dir, outdir_base, natoms, ntypes, ecutwfc, ecutrho, cell_params, atoms, atomic_positions, k_points_grid, metal, degauss=0.0, smearing='gaussian', pressure=0.0):
    """
    Prepare the QE input template with the given parameters.
    """
    template = f"""
&CONTROL
  calculation = '{calculation}'
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
  ecutrho = {ecutrho}"""
    if metal:
        template += f"""
  occupation = 'smearing'
  smearing = '{smearing}'
  degauss = {degauss}"""
    template += f"""
/
&ELECTRONS
  conv_thr = 1.0d-8
  mixing_beta = 0.7
/"""
    if calculation == 'vc-relax' or calculation == 'relax' or calculation == 'md' or calculation == 'vc-md':
       template += f"""
&IONS
  ion_dynamics = 'bfgs'
/
&CELL
  cell_dynamics = 'bfgs'
  press = {pressure}
  press_conv_thr = 0.5  ! kbar
  cell_dofree = 'all'
/"""
       template += f"""
CELL_PARAMETERS {{angstrom}}
{cell_params}
ATOMIC_SPECIES
{atoms}
ATOMIC_POSITIONS {{angstrom}}
{atomic_positions}
K_POINTS {{automatic}}
  {k_points_grid}
"""
    return template

Ge_I = prepare_input(**Ge_I_params, pressure=0.0, calculation='vc-relax')
Ge_II = prepare_input(**Ge_II_params, pressure=100.0, calculation='vc-relax')

if __name__ == "__main__":
    # Print the templates to verify
    print("Ge-I Template:")
    print(Ge_I)
    print("\nGe-II Template:")
    print(Ge_II)
    