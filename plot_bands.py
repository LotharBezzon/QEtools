from urllib import parse
import matplotlib.pyplot as plt
import numpy as np
import argparse

#k_points_path = [r'$\Gamma$', 'X', 'W', 'K', r'$\Gamma$', 'L', 'U', 'W', 'L', 'K|U', '', 'X']
k_points_path = [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'Z', 'P', 'N', r'$Z_1$', 'M|X', '', 'P']

parser = argparse.ArgumentParser(description="Plot band structure from Quantum ESPRESSO bands.dat file.")
parser.add_argument('-i', '--input_prefix', type=str)
parser.add_argument('-p', '--path', action='store_true', required=False,)
parser.add_argument('-n', '--nbnds', type=int, default=12, help='Number of bands to plot (default: 12)')
args = parser.parse_args()

bands_file = f"bands_data/{args.input_prefix}.bands.dat"

with open(bands_file) as f:
    lines = f.readlines()

# Parse header for nbnd and nks
for line in lines:
    if line.strip().startswith('&plot'):
        nbnd = int(line.split('nbnd=')[1].split(',')[0])
        nks = int(line.split('nks=')[1].split('/')[0])
        break

kpoints = []
bands = []

i = 0
while i < len(lines):
    line = lines[i].strip()
    if line and line[0] in "0123456789-." and len(line.split()) == 3:
        # k-point line
        kpt = [float(x) for x in line.split()]
        kpoints.append(kpt)
        # Next line(s): band energies (may be split over two lines)
        band_vals = []
        i += 1
        while len(band_vals) < nbnd:
            vals = [float(x) for x in lines[i].strip().split()]
            band_vals.extend(vals)
            i += 1
        bands.append(band_vals)
    else:
        i += 1

kpoints = np.array(kpoints)        # shape: (nks, 3)
bands = np.array(bands)            # shape: (nks, nbnd)
bands = bands[:, -args.nbnds:]

### Find special k-points
special_kpoints = [(kpoints[0], 0)]

diff_kpoints = np.diff(kpoints, axis=0)
for i in range(1, len(diff_kpoints)):
    if np.linalg.norm(diff_kpoints[i-1] - diff_kpoints[i]) > 1e-4:
        special_kpoints.append((kpoints[i], i))

special_kpoints.append((kpoints[-1], len(kpoints) - 1))

### Define x-axis for plotting
# x starts from 0 and follows the k points path respecting distances in reciprocal space
x = np.linspace(0, np.linalg.norm(special_kpoints[0][0] - special_kpoints[1][0]), special_kpoints[1][1], endpoint=False)
for i in range(1, len(special_kpoints) - 1):
    x = np.concatenate((x, np.linspace(x[-1] + np.diff(x)[-1], x[-1] + np.diff(x)[-1] + np.linalg.norm(special_kpoints[i+1][0] - special_kpoints[i][0]),
                                        special_kpoints[i+1][1] - special_kpoints[i][1], endpoint=False)))
    
x = np.concatenate((x, np.array([x[-1] + np.diff(x)[-1]])))     # shape (nks,)

### Plot the bands
fig = plt.figure(figsize=(8, 5))
plt.plot(x, bands, color='k', linewidth=0.5)
for i in range(len(special_kpoints)):
    plt.axvline(x[special_kpoints[i][1]], color='r', linewidth=0.5, alpha=0.5)

plt.xticks([])
if args.path:
    plt.xticks([x[special_kpoints[i][1]] for i in range(len(special_kpoints))],
               [k_points_path[i] for i in range(len(special_kpoints))])
    
plt.ylabel('Energy (eV)')
    
plt.savefig(f"bands_data/{args.input_prefix}.bands.png", dpi=300, bbox_inches='tight')

plt.show()
