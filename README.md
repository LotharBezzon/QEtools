# QEtools
Some code to work with Quantum Espresso inputs and outputs. It may not be deeply tested.

- **centercell.py**:
Run from the command line with arguments:
  - `-i`: QE input file with atom coordinates to be centered in the cell;
  - `-o` (optional): prefix to be attached to the input file naming the output file (default: 'centered_');
  - `-mc` (optional): save a copy of the input file.


- **templates.py**:
  Contains a function that writes QE input files. Not all calculations may be included (check!). Write a dictionary with desired values as in the examples. Useful to run multiple calculation with different values of some parameters.

- **pressure.py**:
  Run programmatically vc-relaxations at different pressures. The material(s) to be simulated must be described by a dictionary(ies) with its data in `templates.py` (see examples). Pressures at which simulate the material(s) must be specified by the `PRESSURES_KBAR` variable in `pressure.py`. Run with 
  ```bash
  python3 pressure.py -t [prefix material 1] [prefix material 2] ...
  ``` 

- **optimal_parameters.py**:
  Run scf calculations with different values of `ecutwfc` and different k points grids, then plot the resulting energies to check convergence. Need `templates.py` as before. Run with 
  ```bash
  python3 optimal_parameters.py -t [prefix material]
  ``` 
  
- **plot_bands.py**:
  Take as input a `.dat` file from a QE bands calculation. Plot the bands and vertical axes corresponding to the special k points. If you want to print the high-symmetry k points on the x axis use the `-p` option. Run with
  ```bash
  python3 plot_bands.py -i [prefix material] (-p)
  ``` 