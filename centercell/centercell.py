import argparse
import numpy as np
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputfile', help='Input file', type=str)
parser.add_argument('-o', '--outputprefix', help='Output file', default='centered_', type=str, required=False)
parser.add_argument('-mc', '--makecopy', action='store_true', required=False)

args = parser.parse_args()

if args.makecopy==False:
    ans = input('Are you sure you want to overwrite the input file? (y/n) ').lower()
    if ans != 'y':
        sys.exit('If you want to keep the original file, please use the -mc option to keep a copy of the input file')
    
bohr = 0.52917721067    # Bohr to Angstrom conversion factor

with open(args.inputfile, 'r') as f:
    lines = f.readlines()
    for line in lines:
        words = line.split()
        if words[0] == 'ibrav':
            ibrav = int(words[2])
        if words[0] == 'celldm(1)':
            celldm1 = float(words[2])
        if words[0] == 'celldm(2)':
            celldm2 = float(words[2])
        if words[0] == 'celldm(3)':
            celldm3 = float(words[2])
        if words[0] == 'A' and words[1] == '=':
            a = float(words[2])
        if words[0] == 'B' and words[1] == '=':
            b = float(words[2])
        if words[0] == 'C' and words[1] == '=':
            c = float(words[2])
        if words[0] == 'nat':
            nat = int(words[2])

        
        if words[0] == 'ATOMIC_POSITIONS':
            start_coords = lines.index(line)
            

    coords = np.empty((nat, 3))
    for i in range(nat):
        x = float(lines[start_coords + i + 1].split()[1])
        y = float(lines[start_coords + i + 1].split()[2])
        z = float(lines[start_coords + i + 1].split()[3])
        coord = np.array([x, y, z])
        coords[i] = coord

    com = np.mean(coords, axis=0)

    if 'celldm1' in globals():
        a = celldm1 * bohr
    if 'celldm2' in globals():
        b = celldm2 * a
    if 'celldm3' in globals():
        c = celldm3 * a

    if ibrav == 0:      # Free coordinates
        v1 = np.array([a, 0, 0])
        v2 = np.array([0, b, 0])
        v3 = np.array([0, 0, c])
    if ibrav == 1:      # Cubic P
        v1 = np.array([a, 0, 0])
        v2 = np.array([0, a, 0])
        v3 = np.array([0, 0, a])
    if ibrav == 2:      # Cubic F
        v1 = np.array([-a, 0, a]) / 2
        v2 = np.array([0, a, a]) / 2
        v3 = np.array([-a, a, 0]) / 2
    if ibrav == 3:      # Cubic I
        v1 = np.array([a, a, a]) / 2
        v2 = np.array([-a, a, a]) / 2
        v3 = np.array([-a, -a, a]) / 2
    if ibrav == -3:     # Cubic I more symmetric
        v1 = np.array([-a, a, a]) / 2
        v2 = np.array([a, -a, a]) / 2
        v3 = np.array([a, a, -a]) / 2
    if ibrav == 4:      # Hexagonal and trigonal P
        v1 = np.array([a, 0, 0])
        v2 = np.array([-a/2, a*np.sqrt(3)/2, 0])
        v3 = np.array([0, 0, c])


    else:
        raise ValueError('ibrav not supported yet :( \n Please consider contributing to the project at https://github.com/LotharBezzon/QEtools')
    
    cell_center = (v1 + v2 + v3) / 2
    shift = com - cell_center
    coords -= shift

    if args.makecopy == True:
        outfile = args.outputprefix + args.inputfile
    else:
        outfile = args.inputfile

    with open(outfile, 'w') as g:
        for line in lines[:start_coords+1]:
            g.write(line)
        for i in range(nat):
            g.write('C' + '  ' + str(coords[i, 0]) + '  ' + str(coords[i, 1]) + '  ' + str(coords[i, 2]) + '\n')
