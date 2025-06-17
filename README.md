# QEtools
Some code to work with Quantum Espresso inputs and outputs. It may not be deeply tested.
# Centercell
Run from the command line with arguments:
- -i: QE input file with atom coordinates to be centered in the cell;
- -o (optional): prefix to be attached to the input file naming the output file (default: 'centered_');
- -mc (optional): save a copy of the input file.

- **pressure**:
  Run programmatically relaxations at different pressures.